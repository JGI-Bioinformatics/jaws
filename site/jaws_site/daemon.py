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
        self.uploads_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "uploads_subdirectory")
        )
        self.cromwell = Cromwell(conf.get("CROMWELL", "url"))
        self.session = None
        self.active_states = [
            "uploading",
            "upload complete",
            "submitted",
            "queued",
            "running",
            "succeeded",
            "failed",
            "downloading",
        ]
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
        schedule.every(10).seconds.do(self.check_active_runs)
        schedule.every(10).seconds.do(self.update_job_status_logs)
        schedule.every(1).seconds.do(self.send_run_status_logs)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_active_runs(self):
        """
        Check if runs are ready to transition to their next states.
        """
        if self.rpc_client is None:
            self._init_rpc_client()
        self.session = Session()
        try:
            query = (
                self.session.query(Run)
                .filter(Run.status.in_(self.active_states))
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
        for row in query:
            logger.debug(f"Run {row.id}: {row.status}")
            run = Run(row.id, self.session)
            run.check_status()
        self.session.close()

    def update_job_status_logs(self):
        """JTM job status logs are missing some fields; fill them in now."""
        logger.debug("Update job status logs")
        self.session = Session()

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
        self.session.close()

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
            # prepare message
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

            # send via RPC
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

            # mark as sent
            log.sent = True
            try:
                session.commit()
            except Exception as error:
                logger.exception(f"Error updating run_logs as sent: {error}")

        session.close()
