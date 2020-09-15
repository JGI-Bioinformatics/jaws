"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import schedule
import time
import os
import globus_sdk
import logging
from datetime import datetime
from jaws_site.database import Session
from jaws_site.models import Run, Run_Log, Job_Log
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_rpc import rpc_client


logger = logging.getLogger(__package__)


class DataError(Exception):
    pass


class Daemon:
    """
    JAWS Daemon class
    """

    def __init__(self):
        """
        Init obj
        """
        conf = config.conf
        logger.info("Initializing daemon")
        self.site_id = conf.get("SITE", "id")
        self.globus_root_dir = conf.get("GLOBUS", "root_dir")
        self.globus_default_dir = conf.get("GLOBUS", "default_dir")
        self.uploads_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "uploads_subdirectory")
        )
        self.cromwell = Cromwell(conf.get("CROMWELL", "url"))
        self.downloads_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "downloads_subdirectory")
        )
        self.globus_endpoint = conf.get("GLOBUS", "endpoint_id")
        self.downloads_subdir = conf.get("SITE", "downloads_subdirectory")
        self.session = None
        self.operations = {
            "uploading": self.check_if_upload_complete,
            "upload complete": self.submit_run,
            "submitted": self.check_run_cromwell_status,
            "queued": self.check_run_cromwell_status,
            "running": self.check_run_cromwell_status,
            "succeeded": self.transfer_results,
            "failed": self.transfer_results,
            "downloading": self.check_if_download_complete,
        }
        self.terminal_states = ["cancelled", "download complete"]
        self.session = None
        self.rpc_client = None

    def _init_rpc_client(self):
        """Init RPC client for call Central functions"""
        try:
            self.rpc_client = rpc_client.RPC_Client(
                config.conf.get_section("CENTRAL_RPC_CLIENT")
            )
        except Exception as error:
            logger.exception(f"Unable to init central rpc client: {error}")
            raise

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.main_loop)
        schedule.every(1).seconds.do(self.send_run_status_logs)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def _authorize_transfer_client(self, token):
        client_id = config.conf.get("GLOBUS", "client_id")
        client = globus_sdk.NativeAppAuthClient(client_id)
        authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
        return globus_sdk.TransferClient(authorizer=authorizer)

    def main_loop(self):
        """
        Check for runs in particular states.
        """
        if self.rpc_client is None:
            self._init_rpc_client()
        self.session = Session()
        try:
            query = (
                self.session.query(Run)
                .filter(Run.status.in_(self.operations.keys()))
                .all()
            )
        except Exception as error:
            logger.warning(
                f"Failed to select active runs from db: {error}", exc_info=True
            )
            return  # sqlalchemy should reconnect
        num_runs = len(query)
        if num_runs:
            logger.info(f"There are {num_runs} active runs")
        for run in query:
            logger.debug(f"Run {run.id}: {run.status}")
            proc = self.operations.get(run.status, None)
            if proc:
                proc(run)

        # process logs
        self.update_job_status_logs()
        self.session.close()

    def _get_globus_transfer_status(self, run, task_id):
        """
        Query Globus transfer service for transfer task status.
        """
        transfer_rt = run.transfer_refresh_token
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
            task = transfer_client.get_task(task_id)
            globus_status = task["status"]
        except Exception:
            logger.exception("Failed to check Globus upload status", exc_info=True)
            raise
        return globus_status

    def check_if_upload_complete(self, run):
        """
        Check on the status of one uploading run.
        """
        logger.debug(f"Run {run.id}: Check upload status")
        try:
            globus_status = self._get_globus_transfer_status(run, run.upload_task_id)
        except Exception as error:
            logger.exception(f"Failed to check upload {run.upload_task_id}: {error}")
            return
        if globus_status == "FAILED":
            self.update_run_status(run, "upload failed")
        elif globus_status == "INACTIVE":
            self.update_run_status(
                run, "upload inactive", "Your endpoint authorization has expired"
            )
        elif globus_status == "SUCCEEDED":
            self.update_run_status(run, "upload complete")

    def submit_run(self, run):
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {run.id}: Submit to Cromwell")

        # Validate input
        file_path = os.path.join(self.uploads_dir, run.user_id, run.submission_id)
        wdl_file = file_path + ".wdl"
        json_file = file_path + ".json"
        zip_file = file_path + ".zip"  # might not exist
        if not os.path.exists(json_file):
            logger.warning(f"Missing inputs JSON for run {run.id}: {json_file}")
            self.update_run_status(run, "missing input", "Missing inputs JSON file")
            return
        if not os.path.exists(wdl_file):
            logger.warning(f"Missing WDL for run {run.id}: {wdl_file}")
            self.update_run_status(run, "missing input", "Missing WDL file")
            return
        if not os.path.exists(zip_file):
            zip_file = None

        # submit to Cromwell
        cromwell_run_id = self.cromwell.submit(wdl_file, json_file, zip_file)
        if cromwell_run_id:
            run.cromwell_run_id = cromwell_run_id
            self.update_run_status(
                run, "submitted", f"cromwell_run_id={cromwell_run_id}"
            )
        else:
            run.update_run_status(run, "submission failed")

    def check_run_cromwell_status(self, run):
        """
        Check Cromwell for the status of one Run.
        """
        logger.debug(f"Run {run.id}: Check Cromwell status")
        try:
            cromwell_status = self.cromwell.get_status(run.cromwell_run_id)
        except Exception as error:
            logger.error(f"Unable to check Cromwell status of Run {run.id}: {error}")
            return  # try again next time
        logger.debug(f"Run {run.id}: Cromwell status is {cromwell_status}")
        if cromwell_status == "Running":
            if run.status == "submitted":
                try:
                    metadata = self.cromwell.get_metadata(run.cromwell_run_id)
                except Exception:
                    return
                if metadata is None:
                    return
                run.cromwell_workflow_dir = metadata.get("workflowRoot")
                self.update_run_status(run, "queued")
        elif cromwell_status == "Failed":
            self.update_run_status(run, "failed")
        elif cromwell_status == "Succeeded":
            if run.status == "queued":
                self.update_run_status(run, "running")
            if run.status == "running":
                self.update_run_status(run, "succeeded")
        elif cromwell_status == "Aborted":
            if run.status == "queued" or run.status == "running":
                self.update_run_status(run, "cancelled")

    def transfer_results(self, run, metadata=None):
        """
        Send run output via Globus
        """
        logger.debug(f"Run {run.id}: Download output")
        if run.cromwell_workflow_dir is None and metadata is None:
            try:
                metadata = self.cromwell.get_metadata(run.cromwell_run_id)
            except Exception:
                return
            if metadata is None:
                return
            run.cromwell_workflow_dir = metadata.get("workflowRoot")
        transfer_rt = run.transfer_refresh_token
        if not run.cromwell_workflow_dir.startswith(self.globus_root_dir):
            logger.error(
                f"Results dir is not accessible via Globus: {run.cromwell_workflow_dir}"
            )
            return
        rel_cromwell_workflow_dir = os.path.relpath(
            run.cromwell_workflow_dir, self.globus_default_dir
        )
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to get Globus transfer client for run {run.id}", exc_info=True
            )
            return
        try:
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.globus_endpoint,
                run.output_endpoint,
                label=f"Download run {run.id}",
                sync_level="exists",
                verify_checksum=False,
                preserve_timestamp=True,
                notify_on_succeeded=False,
                notify_on_failed=False,
                notify_on_inactive=False,
                skip_activation_check=True,
            )
            if self.globus_root_dir == "/":
                tdata.add_item(
                    run.cromwell_workflow_dir, run.output_dir, recursive=True
                )
            else:
                tdata.add_item(
                    rel_cromwell_workflow_dir, run.output_dir, recursive=True
                )
        except Exception:
            logger.warning(
                f"Failed to prepare download manifest for run {run.id}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to download results with Globus for run {run.id}",
                exc_info=True,
            )
            return
        run.download_task_id = transfer_result["task_id"]
        self.update_run_status(
            run, "downloading", f"download_task_id={run.download_task_id}"
        )

    def check_if_download_complete(self, run):
        """
        If download is complete, change state.
        """
        logger.debug(f"Run {run.id}: Check download status")
        try:
            globus_status = self._get_globus_transfer_status(run, run.download_task_id)
        except Exception as error:
            logger.exception(
                f"Failed to check download {run.download_task_id}: {error}"
            )
            return
        if globus_status == "SUCCEEDED":
            self.update_run_status(run, "download complete")
        elif globus_status == "FAILED":
            self.update_run_status(run, "download failed")

    def update_run_status(self, run, new_status, reason=None):
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        logger.info(f"Run {run.id}: now {new_status}")
        status_from = run.status
        timestamp = datetime.utcnow()

        # update current status in "runs" table
        try:
            run.status = new_status
            run.updated = timestamp
            self.session.commit()
        except Exception as error:
            logger.exception(f"Unable to update Run {run.id}: {error}")

        # populate "reason" field
        if new_status == "submitted":
            reason = f"cromwell_run_id={run.cromwell_run_id}"
        elif new_status == "downloading":
            reason = f"download_task_id={run.download_task_id}"

        # save state transition in "run_logs" table
        try:
            log_entry = Run_Log(
                run_id=run.id,
                status_from=status_from,
                status_to=new_status,
                timestamp=timestamp,
                reason=reason,
            )
            self.session.add(log_entry)
            self.session.commit()
        except Exception as error:
            logger.exception(
                f"Failed to create run_log object for Run {run.id} : {new_status}, {reason}: {error}"
            )
        # notifying Central of state change is handled by send_run_status_logs

    def update_job_status_logs(self):
        """JTM job status logs are missing some fields; fill them in now."""
        logger.debug("Update job status logs")

        # select incomplete job log entries from database
        last_cromwell_run_id = None  # cache last
        last_metadata = None  # no need to get from cromwell repeatedly
        try:
            query = (
                self.session.query(Job_Log)
                .filter_by(task_name=None)
                .filter_by(sent=False)
                .order_by(Job_Log.cromwell_run_id)
            )
        except Exception as error:
            logger.exception(f"Unable to select job_logs: {error}")
            return

        if not query:
            return

        for log in query:
            run_id = log.run_id  # may be None
            cromwell_run_id = log.cromwell_run_id
            cromwell_job_id = log.cromwell_job_id
            logger.debug(f"Job {cromwell_job_id} now {log.status_to}")

            # SELECT run_id FROM RDB
            if not run_id:
                run = (
                    self.session.query(Run)
                    .filter_by(cromwell_run_id=cromwell_run_id)
                    .one_or_none()
                )
                if not run:
                    # JTM may send first job status update before cromwell_run_id is recorded
                    continue
                log.run_id = run_id = run.id

                # update run state if first job to run
                if log.status_to == "running" and run.status == "queued":
                    self.update_run_status(run, "running")

            # TRY TO GET task_name AND attempt FROM CROMWELL METADATA
            if cromwell_run_id == last_cromwell_run_id:
                # another update for same run, reuse previously initialized object
                metadata = last_metadata
            else:
                # get from cromwell
                try:
                    metadata = last_metadata = self.cromwell.get_metadata(
                        cromwell_run_id
                    )
                    last_cromwell_run_id = cromwell_run_id
                except Exception as error:
                    logger.exception(
                        f"Error getting metadata for {cromwell_run_id}: {error}"
                    )
                    continue
            job_info = metadata.get_job_info(cromwell_job_id)
            if job_info is None:
                status = metadata.get("status")
                if status in ("Failed", "Succeeded", "Aborted"):
                    logger.error(
                        f"job_id {cromwell_job_id} not found in inactive run, {cromwell_run_id}"
                    )
                    # If the Run is done and the job_id cannot be found, it's an orphan job
                    # (Cromwell never received the job_id from JTM), so mark as "sent"
                    log.sent = True
                else:
                    # The Cromwell metadata could be a bit outdated; try again next time
                    logger.error(
                        f"job_id {cromwell_job_id} not found in active run, {cromwell_run_id}"
                    )
                continue
            log.attempt = job_info["attempt"]
            log.task_name = job_info["task_name"]
            self.session.commit()

    def send_run_status_logs(self):
        """Send run logs to Central"""

        # get updates from datbase
        try:
            session = Session()
            query = session.query(Run_Log).filter(Run_Log.sent.is_(False)).all()
        except Exception as error:
            logger.exception(f"Unable to select from run_logs: {error}")
            return
        num_logs = len(query)
        if not num_logs:
            return
        logger.debug(f"Sending {num_logs} run logs")

        # send logs via RPC
        for log in query:
            data = {
                "site_id": self.site_id,
                "run_id": log.run_id,
                "status_from": log.status_from,
                "status_to": log.status_to,
                "timestamp": log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                "reason": log.reason,
            }
            # add special fields
            if log.status_to == "submitted":
                run = session.query(Run).get(log.run_id)
                data["cromwell_run_id"] = run.cromwell_run_id
            elif log.status_to == "downloading":
                run = session.query(Run).get(log.run_id)
                data["download_task_id"] = run.download_task_id
            try:
                response = self.rpc_client.request("update_run_logs", data)
            except Exception as error:
                logger.exception(f"RPC update_run_logs error: {error}")
                return
            if "error" in response:
                logger.info(
                    f"RPC update_run_status failed: {response['error']['message']}"
                )
                return
            log.sent = True
            try:
                session.commit()
            except Exception as error:
                logger.exception(f"Error updating run_logs as sent: {error}")
        session.close()
