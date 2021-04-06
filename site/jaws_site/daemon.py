"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import shutil
import schedule
import time
import os
import globus_sdk
import logging
from datetime import datetime
from jaws_site.database import Session
from jaws_site.models import Run, Run_Log
from jaws_site import config, tasks
from jaws_site.cromwell import Cromwell
from jaws_site.globus import GlobusService
from jaws_rpc import rpc_client


logger = logging.getLogger(__package__)


def file_size(file_path: str):
    """
    Checks if a file exists and is a file; if it does, check it's size.

    :param file_path: path to file
    :type file_path: str
    :return: Size in bytes if file exists; else None
    :rtype: int
    """
    if os.path.isfile(file_path):
        return os.path.getsize(file_path)
    else:
        return None


class DataError(Exception):
    pass


class Daemon:
    """
    Daemon that listens to an AMQP messaging queue to coordinate between
    Globus and Cromwell. Daemon will submit to Cromwell upon successful completion of a globus transfer.
    It then tracks the Cromwell status and updates the database depending on the status of the workflow.
    Once completed it will then submit the data back to the submission site.
    """

    def __init__(self):
        logger.info("Initializing daemon")
        self.cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
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
        self.globus = GlobusService()
        self.terminal_states = ["cancelled", "download complete"]
        self.session = None

    def _init_rpc_client(self):
        """
        Instantiates the RPC client to communicate with Central over
        AMQP.
        """
        try:
            central_client = rpc_client.RpcClient(
                config.conf.get_section("CENTRAL_RPC_CLIENT"),
                logger
            )
        except Exception as error:
            logger.exception(f"Unable to init central rpc client: {error}")
            raise
        return central_client

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.main_loop)
        schedule.every(1).seconds.do(self.send_run_status_logs)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def main_loop(self):
        """
        Check for runs in particular states.
        """
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
        self.session.close()

    def _get_globus_transfer_status(self, run, task_id):
        """
        Query Globus transfer service for transfer task status.
        Uploads are done only by JAWS's globus credentials (i.e. confidential app).
        """
        return self.globus.transfer_status(task_id)

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

    def get_uploads_file_path(self, run):
        uploads_dir = config.conf.get("SITE", "uploads_dir")
        return os.path.join(uploads_dir, run.user_id, run.submission_id)

    def get_run_input(self, run):
        """
        Checks if the input files are valid and returns a tuple of their paths.  Otherwise raises an error.

        :param run: The SQLAlchemy ORM model for the run record.
        :type run: obj
        :return: list of [wdl,json,zip] file paths
        :rtype: list
        """
        suffixes_and_required = [("wdl", True), ("json", True), ("zip", False), ("options.json", False)]
        files = []
        file_path = self.get_uploads_file_path(run)
        for (suffix, required) in suffixes_and_required:
            a_file = f"{file_path}.{suffix}"
            a_file_size = file_size(a_file)
            if required:
                if a_file_size is None:
                    raise DataError(f"Input {suffix} file not found: {a_file}")
                elif a_file_size == 0:
                    raise DataError(f"Input {suffix} file is 0-bytes: {a_file}")
            elif a_file_size is not None and a_file_size == 0:
                raise DataError(f"Input {suffix} file is 0-bytes: {a_file}")
            if a_file_size:
                files.append(a_file)
            else:
                files.append(None)
        return files

    def submit_run(self, run):
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {run.id}: Submit to Cromwell")

        try:
            infiles = self.get_run_input(run)
        except Exception as error:
            logger.error(f"Run {run.id}: {error}")
            self.update_run_status(run, "submission failed", f"Bad input: {error}")
            return

        try:
            cromwell_run_id = self.cromwell.submit(*infiles)
        except Exception as error:
            logger.error(f"Run {run.id} submission failed: {error}")
            self.update_run_status(run, "submission failed")
        else:
            run.cromwell_run_id = cromwell_run_id
            self.update_run_status(
                run, "submitted", f"cromwell_run_id={cromwell_run_id}"
            )

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
            if run.status == "submitted":
                self.update_run_status(run, "queued")
            elif run.status == "queued":
                # although Cromwell may consider a Run to be "Running", since it does not distinguish between
                # "queued" and "running", we check the task-log to see if any task is "running"; only once any
                # task is running does the Run transition to the "running" state.
                tasks_status = tasks.get_run_status(self.session, run.id)
                if tasks_status == "running":
                    self.update_run_status(run, "running")
        elif cromwell_status == "Failed":
            self.update_run_status(run, "failed")
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition;
            # sleep is used to avoid multiple state transitions with the same timestamp
            if run.status == "submitted":
                self.update_run_status(run, "queued")
                time.sleep(1)
            if run.status == "queued":
                self.update_run_status(run, "running")
                time.sleep(1)
            if run.status == "running":
                self.update_run_status(run, "succeeded")
        elif cromwell_status == "Aborted":
            self.update_run_status(run, "cancelled")

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
                self.update_run_status(
                    run, "download complete", "No run folder was created"
                )
                return

        file_path = self.get_uploads_file_path(run)
        wdl_file = file_path + ".wdl"
        json_file = file_path + ".json"
        orig_json_file = file_path + ".orig.json"
        zip_file = file_path + ".zip"  # might not exist
        options_file = file_path + ".options.json"  # might not exist
        try:
            shutil.copy(wdl_file, run.cromwell_workflow_dir)
        except Exception as error:
            logger.error(f"Error copying WDL from {wdl_file}->{run.cromwell_workflow_dir}: {error}")
        try:
            shutil.copy(json_file, run.cromwell_workflow_dir)
        except Exception as error:
            logger.error(f"Error copying JSON from {json_file}->{run.cromwell_workflow_dir}: {error}")
        try:
            shutil.copy(orig_json_file, run.cromwell_workflow_dir)
        except Exception as error:
            logger.error(f"Error copying original JSON from {orig_json_file}->{run.cromwell_workflow_dir}: {error}")
        if os.path.exists(zip_file):
            try:
                shutil.copy(zip_file, run.cromwell_workflow_dir)
            except Exception as error:
                logger.error(f"Error copying ZIP from {zip_file}->{run.cromwell_workflow_dir}: {error}")
        if os.path.exists(options_file):
            try:
                shutil.copy(options_file, run.cromwell_workflow_dir)
            except Exception as error:
                logger.error(f"Error copying options-JSON from {options_file}->{run.cromwell_workflow_dir}: {error}")

        transfer_task_id = None
        try:
            transfer_task_id = self.globus.submit_transfer(
                f"Run {run.id}",
                run.output_endpoint,
                run.cromwell_workflow_dir,
                run.output_dir,
            )
        except globus_sdk.GlobusAPIError:
            logger.exception("error while submitting transfer")
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

    def send_run_status_logs(self):
        """Send run logs to Central"""

        # get updates from datbase
        try:
            session = Session()
            query = session.query(Run_Log).filter(Run_Log.sent.is_(False)).all()
        except Exception as error:
            session.close()
            logger.exception(f"Unable to select from run_logs: {error}")
            return
        num_logs = len(query)
        if not num_logs:
            session.close()
            return
        logger.debug(f"Sending {num_logs} run logs")

        central_client = self._init_rpc_client()

        # send logs via RPC
        for log in query:
            data = {
                "site_id": config.conf.get("SITE", "id"),
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
                response = central_client.request("update_run_logs", data)
            except Exception as error:
                logger.exception(f"RPC update_run_logs error: {error}")
                continue
            if "error" in response:
                logger.info(
                    f"RPC update_run_status failed: {response['error']['message']}"
                )
                continue
            log.sent = True
            try:
                session.commit()
            except Exception as error:
                session.rollback()
                logger.exception(f"Error updating run_logs as sent: {error}")
        session.close()
