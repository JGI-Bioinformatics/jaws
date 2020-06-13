"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import schedule
import time
import os
import requests
import globus_sdk
import logging
from datetime import datetime
from jaws_site.database import Session
from jaws_site import models, wfcopy, config


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
        self.workflows_url = conf.get("CROMWELL", "workflows_url")
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
            "cancelling": self.check_run_cromwell_status,
            "succeeded": self.prepare_succeeded_run_output,
            "failed": self.prepare_failed_run_output,
            "ready": self.transfer_results,
            "downloading": self.check_if_download_complete,
        }

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.check_runs)
        while True:
            schedule.run_pending()
            time.sleep(5)

    def _query_by_site(self):
        return (
            self.session.query(models.Run)
            .filter_by(site_id=self.site_id)
            .filter(models.Run.status.in_(self.operations.keys()))
            .all()
        )

    def _query_user_id(self, run):
        return self.session.query(models.User).get(run.user_id)

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
        user = self._query_user_id(run)
        transfer_rt = user.transfer_refresh_token
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
            self.update_run_status(run, "upload inactive")
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
        files = {}
        files["workflowInputs"] = (
            "workflowInputs",
            open(json_file, "r"),
            "application/json",
        )
        files["workflowSource"] = (
            "workflowSource",
            open(wdl_file, "r"),
            "application/json",
        )
        if zip_file:
            files["workflowDependencies"] = (
                "workflowDependencies",
                open(zip_file, "rb"),
                "application/zip",
            )
        try:
            r = requests.post(self.workflows_url, files=files)
        except requests.ConnectionError:
            logger.info("Cromwell server timeout")
            # don't update state; keep trying
            return
        except Exception as err:
            logger.error(f"Unknown error while submitting {run.id} to Cromwell: {err}")
            # don't update state; keep trying
            return
        if r.status_code == 201:
            run.cromwell_id = r.json()["id"]
            self.update_run_status(run, "submitted", "cromwell_id={run.cromwell_id}")
        else:
            logger.warning(f"Run submission {run.id} failed")
            self.update_run_status(run, "submission failed")

    def check_run_cromwell_status(self, run):
        """
        Check Cromwell for the status of one Run.
        """
        try:
            r = requests.get(f"{self.workflows_url}/{run.cromwell_id}/status")
        except requests.ConnectionError:
            logger.warning("Cromwell server timeout")
            return
        if r.status_code != requests.codes.ok:
            return

        cromwell_status = r.json()["status"]
        if cromwell_status == "Running":
            if run.status == "submitted":
                self.update_run_status(run, "running")
        elif cromwell_status == "Failed":
            self.update_run_status(run, "failed")
        elif cromwell_status == "Succeeded":
            self.update_run_status(run, "succeeded")
        elif cromwell_status == "Aborted" and run.status == "cancelling":
            self.update_run_status(run, "cancelled")

    def prepare_failed_run_output(self, run):
        """
        Prepare folder of failed run output.
        """
        try:
            r = requests.get(f"{self.workflows_url}/{run.cromwell_id}/metadata")
        except requests.ConnectionError:
            logger.warning("Cromwell server timeout")
            return
        if r.status_code != requests.codes.ok:
            return
        orig_dir = r.json()["workflowRoot"]
        nice_dir = os.path.join(self.results_dir, str(run.id))
        os.symlink(orig_dir, nice_dir)
        self.update_run_status(run, "ready")

    def prepare_succeeded_run_output(self, run):
        """
        Prepare folder of Run output.
        """
        try:
            r = requests.get(f"{self.workflows_url}/{run.cromwell_id}/metadata")
        except requests.ConnectionError:
            logger.warning("Cromwell server timeout")
            return
        if r.status_code != requests.codes.ok:
            return
        orig_dir = r.json()["workflowRoot"]
        nice_dir = os.path.join(self.results_dir, str(run.id))
        wfcopy.wfcopy(orig_dir, nice_dir)
        self.update_run_status(run, "ready")

    def transfer_results(self, run):
        """
        Send run output via Globus
        """
        user = self._query_user_id(run)
        nice_dir = os.path.join(self.results_dir, str(run.id))
        transfer_rt = user.transfer_refresh_token
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
            log_entry = models.Run_Log(
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
