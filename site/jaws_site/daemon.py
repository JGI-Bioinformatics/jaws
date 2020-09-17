"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import schedule
import time
import os
import globus_sdk
import logging
import re
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
            self.rpc_client = rpc_client.RpcClient(
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
        user = __get_user(run.user_id)
        transfer_rt = user["transfer_refresh_token"]
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
        cromwell_status = None
        if run.cromwell_workflow_dir:
            # all we need is status
            try:
                cromwell_status = self.cromwell.get_status(run.cromwell_run_id)
            except Exception as error:
                logger.error(
                    f"Unable to check Cromwell status of Run {run.id}: {error}"
                )
                return  # try again next time
        else:
            # get complete metadata
            try:
                metadata = self.cromwell.get_metadata(run.cromwell_run_id)
            except Exception:
                return  # try again next time
            if metadata is None:
                return
            run.cromwell_workflow_dir = metadata.get("workflowRoot")
            try:
                self.session.commit()
            except Exception as error:
                self.session.rollback()
                logger.exception(f"Error updating run: {error}")
            cromwell_status = metadata.get("status")

        # check if state has changed; allowed states and transitions for a successful run are:
        # submitted -> queued -> running -> succeeded (with no skipping of states allowed)
        # Additionally, any state can transition directly to cancelled or failed.
        logger.debug(f"Run {run.id}: Cromwell status is {cromwell_status}")
        if cromwell_status == "Running":
            # no skips allowed, so there may be more than one transition
            if run.status == "submitted":
                self.update_run_status(run, "queued")
        elif cromwell_status == "Failed":
            self.update_run_status(run, "failed")
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition
            if run.status == "submitted":
                self.update_run_status(run, "queued")
            if run.status == "queued":
                self.update_run_status(run, "running")
            if run.status == "running":
                self.update_run_status(run, "succeeded")
        elif cromwell_status == "Aborted":
            self.update_run_status(run, "cancelled")

    def __transfer_folder(self, label, transfer_rt, src_dir, dest_endpoint, dest_dir):
        """
        Recursively transfer folder via Globus
        :param label: Label to attach to transfer (e.g. "Run 99")
        :type label: str
        :param transfer_rt: User's Globus transfer refresh token
        :type transfer_rt: str
        :param src_dir: Folder to transfer
        :type src_dir: str
        :param dest_endpoint: Globus endpoint for destination
        :type dest_endpoint: str
        :param dest_dir: Destination path
        :type dest_dir: str
        :return: Globus transfer task id
        :rtype: str
        """
        logger.debug(f"Globus xfer {label}")
        if not src_dir.startswith(self.globus_root_dir):
            logger.error(f"Dir is not accessible via Globus: {src_dir}")
            return None
        rel_src_dir = os.path.relpath(src_dir, self.globus_default_dir)
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to get Globus transfer client to xfer {label}", exc_info=True
            )
            return None
        try:
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.globus_endpoint,
                dest_endpoint,
                label=label,
                sync_level="exists",
                verify_checksum=False,
                preserve_timestamp=True,
                notify_on_succeeded=False,
                notify_on_failed=False,
                notify_on_inactive=False,
                skip_activation_check=True,
            )
            if self.globus_root_dir == "/":
                tdata.add_item(src_dir, dest_dir, recursive=True)
            else:
                tdata.add_item(rel_src_dir, dest_dir, recursive=True)
        except Exception:
            logger.warning(
                f"Failed to prepare download manifest for {label}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to download results with Globus for {label}",
                exc_info=True,
            )
            return None
        return transfer_result["task_id"]

    def transfer_results(self, run, metadata=None):
        """
        Send run output via Globus
        """
        logger.debug(f"Run {run.id}: Download output")
        if run.cromwell_workflow_dir is None:
            if not metadata:
                try:
                    metadata = self.cromwell.get_metadata(run.cromwell_run_id)
                except Exception:
                    return
            if metadata is None:
                return
            outdir = metadata.get("workflowRoot")
            if outdir:
                run.cromwell_workflow_dir = outdir
                try:
                    self.session.commit()
                except Exception as error:
                    self.session.rollback()
                    logger.exception(f"Error updating run: {error}")
            else:
                # This run failed before a folder was created; nothing to xfer
                self.update_run_status(run, "download complete", "No run folder was created")
                return
        transfer_rt = run.transfer_refresh_token
        transfer_task_id = self.__transfer_folder(
            f"Run {run.id}",
            transfer_rt,
            run.cromwell_workflow_dir,
            run.output_endpoint,
            run.output_dir,
        )
        if transfer_task_id:
            run.download_task_id = transfer_task_id
            self.update_run_status(
                run, "downloading", f"download_task_id={run.download_task_id}"
            )
            try:
                self.session.commit()
            except Exception as error:
                self.session.rollback()
                logger.exception(f"Error updating run: {error}")
        # else ignore error; was logged; will try again later

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

    def update_run_status(self, run, new_status, reason=""):
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
            self.session.rollback()
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
            self.session.rollback()
            logger.exception(
                f"Failed to insert log for Run {run.id} : {new_status}, {reason}: {error}"
            )
        # notifying Central of state change is handled by send_run_status_logs

    def update_job_status_logs(self):
        """
        JTM job status logs are missing some fields: run_id, task_name, attempt.
        Fill them in now by querying Cromwell metadata.
        """
        logger.debug("Update job status logs")

        # select incomplete job log entries from database
        last_cromwell_run_id = None  # cache last
        last_metadata = None  # no need to get from cromwell repeatedly
        try:
            query = (
                self.session.query(Job_Log)
                .filter_by(task_name=None)
                .order_by(Job_Log.cromwell_run_id)
            )
        except Exception as error:
            logger.exception(f"Unable to select job_logs: {error}")
            return

        if not query:
            # nothing to do
            return

        for log in query:
            run_id = log.run_id  # may be None
            cromwell_run_id = log.cromwell_run_id
            cromwell_job_id = log.cromwell_job_id
            logger.debug(f"Job {cromwell_job_id} now {log.status_to}")

            # Lookup run_id given cromwell_run_id
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

                try:
                    self.session.commit()
                except Exception as error:
                    self.session.rollback()
                    logger.exception(f"Error updating job log: {error}")

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
                    # (Cromwell never received the job_id from JTM), so set task_name to "ORPHAN"
                    # to mark it as such and so it won't be picked up in the next round.
                    log.task_name = "ORPHAN"
                    try:
                        self.session.commit()
                    except Exception as error:
                        self.session.rollback()
                        logger.exception(f"Error updating job log: {error}")
                    continue
                else:
                    # The Cromwell metadata could be a bit outdated; try again next time
                    logger.debug(
                        f"job_id {cromwell_job_id} not found in active run, {cromwell_run_id}"
                    )
                continue
            log.attempt = job_info["attempt"]
            log.task_name = job_info["task_name"]
            task_dir = job_info[
                "call_root"
            ]  # not saved in db but may be used below for xfer
            try:
                self.session.commit()
            except Exception as error:
                self.session.rollback()
                logger.exception(f"Error updating job log: {error}")

            # if Task complete, then transfer output
            if log.status_to in ["success", "failed"] and task_dir:
                logger.debug(
                    f"Transfer Run {log.run_id}, Task {log.task_name}:{log.attempt}"
                )

                # get Run record
                run = self.session.query(Run).get(log.run_id)

                # set task output dir
                task_output_dir = os.path.normpath(
                    os.path.join(
                        run.output_dir,
                        os.path.relpath(task_dir, run.cromwell_workflow_dir),
                    )
                )

                # set label
                short_task_name = log.task_name.split('.')
                short_task_name = short_task_name[1].split(':')
                label = f"Run {log.run_id} Task {short_task_name[0]}"
                label = re.sub('[^0-9a-zA-Z_]+', ' ', label)

                # recursively transfer task dir
                transfer_task_id = self.__transfer_folder(
                    label,
                    run.transfer_refresh_token,
                    task_dir,
                    run.output_endpoint,
                    task_output_dir,
                )
                log.debug(f"Xfer {log.task_name}: {transfer_task_id}")

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

        if self.rpc_client is None:
            self._init_rpc_client()

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
                session.rollback()
                logger.exception(f"Error updating run_logs as sent: {error}")
        session.close()


def __get_user(user):
    """
    Retrieve a user's information from jaws-user service.
    :param user: The user's user_id
    :type user: str
    :return: User record
    :rtype: dict
    """
    jaws_user_svc = rpc_client.rpc_client(config.conf.get_section("USER_RPC"))
    try:
        response = jaws_user_svc.request("get_user", {"user_id": user})
    except Exception as error:
        logger.exception(f"jaws-user get_user failed: {response['error']['message']}")
    if "error" in response:
        abort(response["error"]["code"], response["error"]["message"])
    return response["result"]
