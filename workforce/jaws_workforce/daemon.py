"""
The workflow_daemon manages the worker pools, resizing as necessary.
"""

import schedule
import time
import os
import logging
from datetime import datetime
from jaws_site.database import Session
from jaws_site.models import Run, Run_Log, Job_Log
from jaws_site import config
from jaws_site.cromwell import Cromwell


logger = logging.getLogger(__package__)


class DataError(Exception):
    pass


class Daemon:
    """
    JAWS Workforce Daemon class
    """

    def __init__(self):
        """
        Init obj
        """
        conf = config.conf
        logger.info("Initializing daemon")
        self.site_id = conf.get("SITE", "id")
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

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(1).minutes.do(self.eval_workforce)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def eval_workforce(self):
        """
        Compare workforce to queue and reize as appropriate.
        """
        session = Session()
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

        # resize
        self.resize_workforce()
        session.close()


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
            query = self.session.query(Run_Log).filter(Run_Log.sent.is_(False)).all()
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
                run = self.session.query(Run).get(log.run_id)
                data["cromwell_run_id"] = run.cromwell_run_id
            elif log.status_to == "downloading":
                run = self.session.query(Run).get(log.run_id)
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
                self.session.commit()
            except Exception as error:
                logger.exception(f"Error updating run_logs as sent: {error}")
