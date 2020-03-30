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

from jaws_site import models, wfcopy
from jaws_site.config import jaws_config


class DataError(Exception):
    pass


class JAWSd:
    """
    JAWS Daemon class
    """

    interval_seconds = 15
    watched_states = (
        "uploading",
        "staged",
        "succeeded",
        "ready",
        "downloading",
        "submitted",
        "queued",
        "running",
        "aborting",
    )

    def __init__(self, db):
        """
        Init obj
        """
        self.conf = jaws_config
        self.logger = logging.getLogger(__package__)
        self.logger.debug("Initializing daemon")
        self.site_id = self.conf.get("SITE", "id")
        self.session = db.session()
        self.staging_dir = os.path.join(
            self.conf.get("GLOBUS", "root_dir"),
            self.conf.get("SITE", "staging_subdirectory")
        )
        self.workflows_url = self.conf.get("CROMWELL", "workflows_url")
        self.results_dir = os.path.join(
            self.conf.get("GLOBUS", "root_dir"),
            self.conf.get("SITE", "results_subdirectory")
        )
        self.globus_endpoint = self.conf.get("GLOBUS", "endpoint_id")
        self.results_subdir = self.conf.get("SITE", "results_subdirectory")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(self.interval_seconds).seconds.do(self.check_runs)
        while True:
            schedule.run_pending()
            time.sleep(5)

    def _query_by_site(self):
        return (
            self.session.query(models.Run)
            .filter_by(site_id=self.site_id)
            .filter(models.Run.status.in_(self.watched_states))
            .all()
        )

    def _query_user_id(self, run):
        return self.session.query(models.User).get(run.user_id)

    def _authorize_transfer_client(self, token):
        client_id = self.conf.get("GLOBUS", "client_id")
        client = globus_sdk.NativeAppAuthClient(client_id)
        authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
        return globus_sdk.TransferClient(authorizer=authorizer)

    def check_runs(self):
        """
        Check for runs in particular states.
        """
        self.logger.debug("Query db for updated runs")
        try:
            q = self._query_by_site()
        except Exception:
            self.logger.warning(
                "Failed to query db for runs in watched states", exc_info=True
            )
            return  # sqlalchemy should reconnect
        for run in q:
            if run.status == "uploading":
                self.check_transfer_status(run)
            elif run.status == "staged":
                self.submit_run(run)
            elif run.status == "succeeded":
                self.prepare_run_output(run)
            elif run.status == "ready":
                self.transfer_results(run)
            elif run.status == "downloading":
                self.check_if_download_complete(run)
            else:
                self.check_run_status(run)
        self.session.close_all()

    def check_transfer_status(self, run):
        """
        Check on the status of one uploading run.
        """
        user = self._query_user_id(run)
        transfer_rt = user.transfer_refresh_token
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
            task = transfer_client.get_task(run.upload_task_id)
            globus_status = task["status"]
        except Exception:
            self.logger.warning("Failed to check Globus upload status", exc_info=True)
            raise
        if globus_status == "FAILED":
            run.status = "upload failed"
            self.session.commit()
        elif globus_status == "INACTIVE":
            run.status = "upload inactive"
            self.session.commit()
        elif globus_status == "SUCCEEDED":
            self.submit_run(run)

    def submit_run(self, run):
        """
        Submit a run to Cromwell.
        """
        file_path = os.path.join(self.staging_dir, run.user_id,
                                 run.submission_id)
        wdl_file = file_path + ".wdl"
        json_file = file_path + ".json"
        zip_file = file_path + ".zip"  # might not exist
        files = {}
        try:
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
            if os.path.exists(zip_file):
                files["workflowDependencies"] = (
                    "workflowDependencies",
                    open(zip_file, "rb"),
                    "application/zip",
                )
        except Exception:
            self.logger.warning(f"Invalid input for run {run.id}", exc_info=True)
            run.status = "invalid_input"
            self.session.commit()
        try:
            r = requests.post(self.workflows_url, files=files)
        except requests.ConnectionError:
            self.logger.info("Cromwell server timeout")
            # don't update state; keep trying
            return
        if r.status_code == 201:
            run.cromwell_id = r.json()["id"]
            run.status = "submitted"
        else:
            self.logger.warning(f"Run submission {run.id} failed")
            run.status = "submission failed"
        self.session.commit()

    def check_run_status(self, run):
        """
        Check on the status of one Run.
        """
        try:
            r = requests.get("%s/%s/status" % (self.workflows_url,
                                               run.cromwell_id))
        except requests.ConnectionError:
            self.logger.warning("Cromwell server timeout")
            return
        if r.status_code != requests.codes.ok:
            return

        cromwell_status = r.json()["status"]
        if cromwell_status == "Running":
            if run.status == "submitted":
                run.status = "running"
                self.session.commit()
        elif cromwell_status == "Failed":
            run.status = "failed"
            self.session.commit()
        elif cromwell_status == "Succeeded":
            run.status = "succeeded"
            self.session.commit()
            self.prepare_run_output(run)
        elif cromwell_status == "Aborted" and run.status == "aborting":
            run.status = "aborted"
            self.session.commit()

    def prepare_run_output(self, run):
        """
        Prepare folder of Run output.
        """
        url = "%s/%s/metadata" % (self.workflows_url, run.cromwell_id)
        try:
            r = requests.get(url)
        except requests.ConnectionError:
            self.logger.warning("Cromwell server timeout")
            return
        if r.status_code != requests.codes.ok:
            return
        orig_dir = r.json()["workflowRoot"]
        nice_dir = os.path.join(self.results_dir, str(run.id))
        wfcopy.wfcopy(orig_dir, nice_dir, flattenShardDir=False, verbose=False)
        run.status = "ready"
        self.session.commit()
        self.transfer_results(run)

    def transfer_results(self, run):
        """
        Send run output via Globus
        """
        user = self._query_user_id(run)
        nice_dir = os.path.join(self.results_subdir, str(run.id))
        transfer_rt = user.transfer_refresh_token
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
        except globus_sdk.GlobusAPIError:
            self.logger.warning(
                f"Failed to get Globus transfer client for run {run.id}", exc_info=True
            )
            return
        try:
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.globus_endpoint,
                run.output_endpoint,
                label=f"run_id={run.id}",
                sync_level="checksum",
                verify_checksum=True,
                preserve_timestamp=True,
                notify_on_succeeded=True,
                notify_on_failed=True,
                notify_on_inactive=True,
                skip_activation_check=True,
            )
            tdata.add_item(nice_dir, run.output_dir, recursive=True)
        except Exception:
            self.logger.warning(
                f"Failed to prepare download manifest for run {run.id}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except globus_sdk.GlobusAPIError:
            self.logger.warning(
                f"Failed to download results with Globus for run {run.id}",
                exc_info=True,
            )
            return
        run.download_task_id = transfer_result["task_id"]
        self.session.commit()

    def check_if_download_complete(self, run):
        """
        If download is complete, change state.
        """
        user = self._query_user_id(run)
        task_id = run.download_task_id
        token = user.transfer_refresh_token
        try:
            transfer_client = self._authorize_transfer_client(token)
            result = transfer_client.get_task(task_id)
            status = result["status"]
        except globus_sdk.GlobusAPIError:
            self.logger.warning(
                f"Failed to get Globus Transfer Client for run {run.id}"
            )
            raise
        if status == "SUCCEEDED":
            run.status = "finished"
            self.session.commit()
        elif status == "FAILED":
            run.status = "download failed"
            self.session.commit()
