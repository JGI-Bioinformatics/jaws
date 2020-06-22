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
from jaws_site import wfcopy, config
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
        self.staging_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "staging_subdirectory")
        )
        self.cromwell = Cromwell(conf.get("CROMWELL", "url"))
        self.results_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "results_subdirectory")
        )
        self.globus_endpoint = conf.get("GLOBUS", "endpoint_id")
        self.results_subdir = conf.get("SITE", "results_subdirectory")
        self.session = None
        self.operations = {
            "uploading": self.check_if_upload_complete,
            "upload succeeded": self.submit_run,
            "staged": self.submit_run,
            "submitted": self.check_run_cromwell_status,
            "queued": self.check_run_cromwell_status,
            "running": self.check_run_cromwell_status,
            "succeeded": self.prepare_run_output,
            "failed": self.prepare_run_output,
            "ready": self.transfer_results,
            "downloading": self.check_if_download_complete,
        }
        self.terminal_states = ["cancelled", "finished"]
        self.session = None
        self.rpc_client = None

    def init_rpc_client(self):
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
        while True:
            schedule.run_pending()
            time.sleep(5)

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
            logger.warning(f"Failed to select runs from db: {error}", exc_info=True)
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
        self.send_run_status_logs()
        self.update_job_status_logs()
        self.send_job_status_logs()
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
            self.update_run_status(run, "upload succeeded")

    def submit_run(self, run):
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {run.id}: Submit to Cromwell")

        # Validate input
        file_path = os.path.join(self.staging_dir, run.user_id, run.submission_id)
        wdl_file = file_path + ".wdl"
        json_file = file_path + ".json"
        zip_file = file_path + ".zip"  # might not exist
        if not os.path.exists(json_file):
            logger.warning(f"Missing inputs JSON for run {run.id}: {json_file}")
            self.update_run_status(run, "invalid input", "Missing inputs JSON file")
            return
        if not os.path.exists(wdl_file):
            logger.warning(f"Missing WDL for run {run.id}: {wdl_file}")
            self.update_run_status(run, "invalid input", "Missing WDL file")
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
                self.update_run_status(run, "queued")
        elif cromwell_status == "Failed":
            if run.status == "queued":
                self.update_run_status(run, "running")
            if run.status == "running":
                self.update_run_status(run, "failed")
        elif cromwell_status == "Succeeded":
            if run.status == "queued":
                self.update_run_status(run, "running")
            if run.status == "running":
                self.update_run_status(run, "succeeded")
        elif cromwell_status == "Aborted":
            if run.status == "queued" or run.status == "running":
                self.update_run_status(run, "cancelled")

    def prepare_run_output(self, run, metadata=None):
        """
        Prepare folder of Run output.
        """
        logger.debug(f"Run {run.id}: Prepare output")
        if metadata is None:
            try:
                metadata = self.cromwell.get_metadata(run.cromwell_run_id)
            except Exception:
                return
            if metadata is None:
                return
        logger.info(f"Run {run.id}: Prepare output of successful run")
        orig_dir = metadata.get("workflowRoot")
        nice_dir = os.path.join(self.results_dir, str(run.id))
        if run.status == "succeeded":
            wfcopy.wfcopy(orig_dir, nice_dir)
        else:  # failed
            os.symlink(orig_dir, nice_dir)
        self.update_run_status(run, "ready")

    def transfer_results(self, run):
        """
        Send run output via Globus
        """
        logger.debug(f"Run {run.id}: Download output")
        nice_dir = os.path.join(self.results_dir, str(run.id))
        transfer_rt = run.transfer_refresh_token
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
                sync_level="checksum",
                verify_checksum=True,
                preserve_timestamp=True,
                notify_on_succeeded=False,
                notify_on_failed=True,
                notify_on_inactive=True,
                skip_activation_check=True,
            )
            tdata.add_item(nice_dir, run.output_dir, recursive=True)
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
            self.update_run_status(run, "finished")
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

        # record log state transition in "run_logs"
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

        # notify Central
        log_msg = {
            "run_id": run.id,
            "status": new_status,
            "timestamp": timestamp.strftime("%Y-%m-%d %H:%M:%S"),
        }
        if new_status == "submitted":
            log_msg["cromwell_run_id"] = run.cromwell_run_id
        elif new_status == "downloading":
            log_msg["download_task_id"] = run.download_task_id
        response = self.rpc_client.request("update_run_status", log_msg)
        if "error" in response:
            logger.info(f"RPC update_run_status failed: {response['error']['message']}")
            # failure to notify Central doesn't block state transition
        elif "result" in response:
            logger.debug("Central notified of change of run state")

        # delete finished runs from db
        if new_status in self.terminal_states:
            Run.query.filter_by(id=run.id).delete()
            self.session.commit()

    def update_job_status_logs(self):
        """JTM job status logs are missing some fields; fill them in now."""
        logger.debug("Update job status logs")

        # select incomplete job log entries from database
        last_cromwell_run_id = None  # cache last
        last_metadata = None  # no need to get from cromwell repeatedly
        try:
            query = self.session.query(Job_Log).filter_by(task_name=None)
        except Exception as error:
            logger.exception(f"Unable to select job_logs: {error}")
            return
        for log in query:
            run_id = log.run_id
            cromwell_run_id = log.cromwell_run_id
            cromwell_job_id = log.cromwell_job_id

            # SELECT run_id FROM RDB
            run = (
                self.session.query(Run)
                .filter_by(cromwell_run_id=cromwell_run_id)
                .one_or_none()
            )
            if not run:
                # JTM may send first job status update before cromwell_run_id is recorded
                logger.debug(f"cromwell workflow run not found: {cromwell_run_id}")
                continue

            # update run state if first job to run
            if log.status_to == "running" and run.status == "queued":
                self.update_run_status(run, "running")

            # update log entry with run_id
            logger.debug(f"Cromwell run {cromwell_job_id} => JAWS run {run.id}")
            log.run_id = run.id

            # TRY TO GET task_name AND attempt FROM CROMWELL METADATA
            logger.debug(f"Run {run_id}: Get job info for job {cromwell_job_id}")
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
                logger.error(
                    f"job_id {cromwell_job_id} not found in {cromwell_run_id}: {job_info}"
                )
                continue
            log.attempt = job_info["attempt"]
            log.task_name = job_info["task_name"]
            self.session.commit()

    def send_run_status_logs(self):
        """Batch and send run logs to Central"""
        logger.debug("Send run log updates to Central")

        # get updates from datbase
        try:
            query = self.session.query(Run_Log)
        except Exception as error:
            logger.exception(f"Unable to select from run_logs: {error}")
            return

        logs = []
        for log in query:
            log_entry = [
                log.run_id,
                log.status_from,
                log.status_to,
                log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                log.reason,
            ]
            logs.append(log_entry)

        num_rows = len(logs)
        if num_rows == 0:
            return

        # notify Central
        logger.debug(f"RPC update_run_logs: {num_rows} items")
        data = {"logs": logs, "site_id": self.site_id}
        try:
            response = self.rpc_client.request("update_run_logs", data)
        except Exception as error:
            logger.exception(f"RPC update_run_logs error: {error}")
            return
        if "error" in response:
            logger.info(f"RPC update_run_status failed: {response['error']['message']}")
            return
        try:
            query.delete(synchronize_session=False)
            self.session.commit()
        except Exception as error:
            logger.exception(f"Error deleting rows from run_logs: {error}")

    def send_job_status_logs(self):
        """Batch and send job logs to Central"""
        logger.info("Checking for new job status log entries")

        # get updates from datbase
        try:
            query = self.session.query(Job_Log).filter(Job_Log.run_id.isnot(None))
        except Exception as error:
            logger.error(f"Unable to select job_logs: {error}")
            return

        # prepare message
        logs = []
        for log in query:
            run_id = log.run_id
            cromwell_run_id = log.cromwell_run_id
            cromwell_job_id = log.cromwell_job_id

            # get task_name and attempt from cromwell metadata
            logger.debug(f"Run {run_id}: Get job info for job {cromwell_job_id}")
            try:
                metadata = self.cromwell.get_metadata(cromwell_run_id)
            except Exception as error:
                logger.error(f"Error getting metadata for {cromwell_run_id}: {error}")
                continue
            job_info = metadata.get_job_info(cromwell_job_id)
            if job_info is None:
                logger.error(
                    f"job_id {cromwell_job_id} not found in {cromwell_run_id}: {job_info}"
                )
                continue

            # append new log entry
            log_entry = [
                log.run_id,
                job_info["task_name"],
                int(job_info["attempt"]),
                log.cromwell_job_id,
                log.status_from,
                log.status_to,
                log.cromwell_job_id,
                log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                log.reason,
            ]
            logs.append(log_entry)

        num_rows = len(logs)
        if num_rows == 0:
            return

        # send message to Central
        logger.info(f"Sending {num_rows} job status updates to Central")
        data = {"logs": logs, "site": self.site_id}
        try:
            response = self.rpc_client.request("update_job_logs", data)
        except Exception as error:
            logger.exception(f"RPC update_job_logs error: {error}")
            return
        if "error" in response:
            logger.error(f"RPC update_job_logs failed: {response['error']['message']}")
            return
        # delete logs from database
        try:
            query.delete(synchronize_session=False)
            self.session.commit()
        except Exception as error:
            logger.error(f"Failed to delete job logs: {error}")
