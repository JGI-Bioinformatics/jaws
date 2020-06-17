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


logger = None  # set by Daemon init


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
        global logger
        logger = logging.getLogger(__package__)
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
        self.rpc_client = None

    def init_rpc_client(self):
        """Init RPC client for call Central functions"""
        try:
            self.rpc_client = rpc_client.RPC_Client(
                config.conf.get_section("CENTRAL_RPC_CLIENT")
            )
        except Exception as error:
            logger.error(f"Unable to init central rpc client: {error}")
            raise

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.check_runs)
        schedule.every(10).seconds.do(self.process_logs)
        while True:
            schedule.run_pending()
            time.sleep(5)

    def _query_by_site(self):
        return (
            self.session.query(Run)
            .filter_by(site_id=self.site_id)
            .filter(Run.status.in_(self.operations.keys()))
            .all()
        )

    def _authorize_transfer_client(self, token):
        client_id = config.conf.get("GLOBUS", "client_id")
        client = globus_sdk.NativeAppAuthClient(client_id)
        authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
        return globus_sdk.TransferClient(authorizer=authorizer)

    def check_runs(self):
        """
        Check for runs in particular states.
        """
        self.session = Session()
        try:
            query = self._query_by_site()
        except Exception:
            logger.warning(
                "Failed to query db for runs in watched states", exc_info=True
            )
            return  # sqlalchemy should reconnect
        num_runs = len(query)
        if num_runs:
            logger.debug(f"Checking on {num_runs} active runs")
        for run in query:
            proc = self.operations.get(run.status, None)
            if proc:
                proc(run)
        self.session.close_all()

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
        try:
            globus_status = self._get_globus_transfer_status(run, run.upload_task_id)
        except Exception:
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
        logger.info(f"Run {run.id}: Submit to Cromwell")

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
        try:
            cromwell_status = self.cromwell.get_status(run.cromwell_run_id)
        except Exception:
            return  # try again next time
        if cromwell_status == "Running":
            if run.status == "submitted":
                self.update_run_status(run, "running")
        elif cromwell_status == "Failed":
            if run.status == "running":
                self.update_run_status(run, "failed")
        elif cromwell_status == "Succeeded":
            if run.status == "running":
                self.update_run_status(run, "succeeded")
        elif cromwell_status == "Aborted":
            if run.status == "running":
                self.update_run_status(run, "cancelled")

    def prepare_run_output(self, run, metadata=None):
        """
        Prepare folder of Run output.
        """
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
        try:
            globus_status = self._get_globus_transfer_status(run, run.download_task_id)
        except Exception:
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
        run.status = new_status
        run.updated = timestamp
        self.session.commit()

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
            "timestamp": timestamp.strftime("%Y-%m-%d %H:%M:%S")
        }
        response = self.rpc_client.request("update_run_status", log_msg)
        if "error" in response:
            logger.info(f"RPC update_run_status failed: {response['error']['message']}")
            # failure to notify Central doesn't block state transition

    def process_logs(self):
        """Periodically process logs"""
        self.send_run_status_logs()
        self.update_job_status_logs()
        self.send_job_status_logs()

    def update_job_status_logs(self):
        """JTM job status logs are missing some fields; fill them in now."""

        # select incomplete job log entries from database
        last_cromwell_run_id = (
            None  # there are often many job log entries for the same run;
        )
        last_metadata = None  # no need to get from cromwell repeatedly
        self.session = Session()
        query = (
            self.session.query(Job_Log)
            .filter_by(task_name=None)
            .order_by(Job_Log.cromwell_id)
        )
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
                logger.error(f"cromwell workflow run not found: {cromwell_run_id}")
                continue

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
                    logger.error(
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
        if not self.rpc_client:
            return
        logger.debug("Checking for new run status logs")

        # get updates from datbase
        logs = []
        self.session = Session()
        query = self.session.query(Run_Log)
        for log in query:
            timestamp = log.timestamp.strftime("%Y-%m-%d %H:%M:%S")
            log_entry = {
                "run_id": log.run_id,
                "status_from": log.status_from,
                "status_to": log.status_to,
                "timestamp": timestamp,
                "reason": log.reason,
            }
            logs.append(log_entry)

        # notify Central
        response = self.rpc_client.request("update_run_logs", logs)
        if "error" in response:
            logger.info(f"RPC update_run_status failed: {response['error']['message']}")
            self.session.rollback()
        else:
            query.delete(synchronize_session=False)
            self.session.commit()

    def send_job_status_logs(self):
        """Batch and send job logs to Central"""
        if not self.rpc_client:
            return
        logger.info("Checking for new job status log entries")

        # get updates from datbase
        logs = []
        self.session = Session()
        query = self.session.query(Job_Log)
        for log in query:
            run_id = log.run_id
            cromwell_run_id = log.cromwell_run_id
            cromwell_job_id = log.cromwell_job_id
            timestamp = log.timestamp.strftime("%Y-%m-%d %H:%M:%S")
            item = {
                "run_id": log.run_id,
                "cromwell_job_id": log.cromwell_job_id,
                "status_from": log.status_from,
                "status_to": log.status_to,
                "timestamp": timestamp,
                "reason": log.reason,
            }

            # TRY TO GET task_name AND attempt FROM CROMWELL METADATA
            logger.debug(f"Run {run_id}: Get job info for job {cromwell_job_id}")
            try:
                metadata = self.cromwell.get_metadata(cromwell_run_id)
            except Exception as error:
                logger.error(f"Error getting metadata for {cromwell_run_id}: {error}")
                continue
            job_info = metadata.get_job_info(cromwell_job_id)
            if job_info is None:
                logger.error(f"job_id {cromwell_job_id} not found in {cromwell_run_id}: {job_info}")
                continue
            item["attempt"] = job_info["attempt"]
            item["task_name"] = job_info["task_name"]

            # append new log entry
            logs.append(item)
        num_rows = len(logs)
        if num_rows == 0:
            return

        # send message to Central
        data = {"logs": logs, "site": self.site_id}
        logger.info(f"Sending {num_rows} job status updates to Central")
        result = self.rpc_client.request("update_job_logs", data)
        if "result" in result:
            # delete logs from database
            query.delete(synchronize_session=False)
            self.session.commit()
