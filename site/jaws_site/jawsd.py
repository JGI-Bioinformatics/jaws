"""
JAWS Daemon process periodically checks on runs and performs actions to usher them to the next state.
"""

import schedule
import time
import os
import requests
import globus_sdk
from globus_sdk import TransferClient, RefreshTokenAuthorizer, GlobusAPIError
import logging
import click

from jaws_site import models, wfcopy, config, database


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
    site_id = None
    session = None
    staging_dir = None
    workflows_url = None
    results_dir = None
    globus_endpoint = None
    results_subdir = None

    def __init__(self, conf, db):
        """
        Init obj
        """
        self.logger = logging.getLogger(__package__)
        self.logger.debug("Initializing daemon")
        self.site_id = conf.get("SITE", "id")
        self.session = db.session()
        self.staging_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "staging_subdirectory")
        )
        self.workflows_url = conf.get("CROMWELL", "workflows_url")
        self.results_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "results_subdirectory")
        )
        self.globus_endpoint = conf.get("GLOBUS", "endpoint_id")
        self.results_subdir = conf.get("SITE", "results_subdirectory")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(self.interval_seconds).seconds.do(self.check_runs)
        while True:
            schedule.run_pending()
            time.sleep(5)

    def check_runs(self):
        """
        Check for runs in particular states.
        """
        self.logger.debug("Query db for updated runs")
        try:
            q = (
                self.session.query(models.Run)
                .filter_by(site_id=self.site_id)
                .filter(models.Run.status.in_(self.watched_states))
                .all()
            )
        except Exception:
            self.logger.warning(
                "Failed to query db for runs in watched states", exc_info=True
            )
            return  # sqlalchemy should reconnect and the connection shall be available next time
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
        user = self.session.query(models.User).get(run.user_id)
        try:
            authorizer = RefreshTokenAuthorizer(user.transfer_refresh_token)
            transfer_client = TransferClient(authorizer)
            task = transfer_client.get_task(run.upload_task_id)
            globus_status = task["status"]
        except Exception:
            self.logger.warning("Failed to check Globus upload status", exc_info=True)
            return
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
        wdl_file = os.path.join(
            self.staging_dir, run.user_id, run.submission_id + ".wdl"
        )
        json_file = os.path.join(
            self.staging_dir, run.user_id, run.submission_id + ".json"
        )
        zip_file = os.path.join(
            self.staging_dir, run.user_id, run.submission_id + ".zip"
        )  # might not exist
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
            run.status = "invalid input"
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
            self.logger.warn(f"Run submission {run.id} failed")
            run.status = "submission failed"
        self.session.commit()

    def check_run_status(self, run):
        """
        Check on the status of one Run.
        """
        try:
            r = requests.get("%s/%s/status" % (self.workflows_url, run.cromwell_id))
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
        user = self.session.query(models.User).get(run.user_id)
        nice_dir = os.path.join(self.results_subdir, str(run.id))
        token = user.transfer_refresh_token
        try:
            authorizer = RefreshTokenAuthorizer(token)
            transfer_client = TransferClient(authorizer=authorizer)
        except Exception:
            self.logger.warning(
                f"Failed to get Globus transfer client for run {run.id}", exc_info=True
            )
            return
        label = "run_id=%s" % (run.id,)
        try:
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.globus_endpoint,
                run.dest_endpoint,
                label=label,
                sync_level="checksum",
                verify_checksum=True,
                preserve_timestamp=True,
                notify_on_succeeded=True,
                notify_on_failed=True,
                notify_on_inactive=True,
                skip_activation_check=True,
            )
            tdata.add_item(nice_dir, run.dest_path, recursive=True)
        except Exception:
            self.logger.warning(
                f"Failed to prepare download manifest for run {run.id}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except GlobusAPIError:
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
        user = self.session.query(models.User).get(run.user_id)
        task_id = run.download_task_id
        user = self.session.query(models.User).get(run.user_id)
        token = user.transfer_refresh_token
        authorizer = RefreshTokenAuthorizer(token)
        try:
            transfer_client = TransferClient(authorizer=authorizer)
        except Exception:
            self.logger.warning(
                f"Failed to get Globus Transfer Client for run {run.id}"
            )
        try:
            result = transfer_client.get_task(task_id)
            status = result["status"]
        except Exception:
            self.logger.warning(f"Failed to get download status for task {task_id}")
        if status == "SUCCEEDED":
            run.status = "finished"
            self.session.commit()
        elif status == "FAILED":
            run.status == "download failed"
            self.session.commit()


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """JAWS-Site RPC functions"""
    pass


@cli.command()
def check_runs():
    """Check if any runs need action."""
    conf = config.conf
    db = database.db
    jawsd = JAWSd(conf, db)
    jawsd.check_runs()


@cli.command()
@click.argument("user_id")
@click.argument("submission_id")
@click.argument("dest_endpoint")
@click.argument("dest_dir")
def submit_run(user_id: str, submission_id: str, dest_endpoint: str, dest_dir: str):
    """Submit a new run."""
    conf = config.conf
    db = database.db
    jawsd = JAWSd(conf, db)
    user = jawsd.session.query(models.User).get(user_id).one_or_none()
    if user is None:
        raise DataError(f"User does not exist: {user_id}")
    run = models.Run(
        user_id=user_id,
        status="staged",
        submission_id=submission_id,
        upload_task_id="NA",
        globus_endpoint=dest_endpoint,
        outdir=dest_dir
    )
    jawsd.session.add(run)
    jawsd.session.commit()
    jawsd.submit_run(run)
