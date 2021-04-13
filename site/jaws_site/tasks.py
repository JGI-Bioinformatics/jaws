"""
This class represents the task-logs which are stored in a rdb and referenced by cromwell_run_id instead of jaws run_id.
"""

import logging
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from jaws_site.models import Job_Log, Run
from jaws_site import config
from jaws_site.cromwell import Cromwell


logger = logging.getLogger(__package__)


job_status_value = {
    "created": 0,
    "ready": 1,
    "queued": 2,
    "pending": 3,
    "running": 4,
    "success": 5,
    "failed": 6,
}


class TaskLogDbError(Exception):
    """This is raised on database errors."""

    pass


class TaskLogRunNotFoundError(Exception):
    """This is raised when the requested JAWS run_id is not found in the rdb."""

    pass


class TaskLogError(Exception):
    pass


class TaskLog:
    """
    The task log is a record of state transitions for tasks.
    """

    def __init__(self, session):
        self.cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
        self.session = session

    def _save_job_log(self, job_log):
        """
        Insert or ignore job_log object into RDb.
        This ignores redundant entries because JTM may send duplicate messages.
        :param job_log: sqlalchemy model for the job log
        :type job_log: jaws_site.models.Job_Log
        :return: table of job logs
        :rtype: list
        """
        try:
            self.session.add(job_log)
            self.session.commit()
            logger.debug(f"Job {job_log.cromwell_job_id} status saved")
        except IntegrityError:
            # JTM sometimes sends duplicate messages; ignore
            self.session.rollback()
            logger.warning(
                f"Job {job_log.cromwell_job_id} status is duplicate; ignored"
            )
        except Exception as error:
            self.session.rollback()
            raise TaskLogDbError(
                f"Error saving job log, {job_log.cromwell_job_id}:{job_log.status_to}, in db: {error}"
            )

    def save_job_log(
        self,
        cromwell_run_id: str,
        cromwell_job_id: int,
        status_from: str,
        status_to: str,
        timestamp: str,
        reason: str = None,
    ):
        """Save a task state transition log"""
        logger.debug(
            f"Save job log: {cromwell_run_id}:{cromwell_job_id}:{status_from}:{status_to}:{reason}"
        )
        try:
            job_log = Job_Log(
                cromwell_run_id=cromwell_run_id,
                cromwell_job_id=cromwell_job_id,
                status_from=status_from,
                status_to=status_to,
                timestamp=timestamp,
                reason=reason,
            )
        except Exception as error:
            raise TaskLogError(
                f"Error initializing job_log {cromwell_job_id}:{status_to}: {error}"
            )
        self._save_job_log(job_log)

    def _get_job_logs(self, cromwell_run_ids):
        """Get all logs associated with the query cromwell run id.
        :param cromwell_run_ids: List of cromwell run ids of main workflow and any subworkflows
        :type cromwell_run_ids: list
        :return: table of job logs
        :rtype: list
        """
        try:
            table = (
                self.session.query(Job_Log)
                .filter(Job_Log.cromwell_run_id.in_(cromwell_run_ids))
                .all()
            )
        except SQLAlchemyError as error:
            raise TaskLogError(
                f"SQLAlchemy error retrieving task logs from db for run {cromwell_run_ids}: {error}"
            )
        except Exception as error:
            raise TaskLogError(
                f"Unknown error retrieving task logs from db for run {cromwell_run_ids}: {error}"
            )
        result = []
        for row in table:
            log_entry = [
                row.cromwell_run_id,
                row.cromwell_job_id,
                row.status_from,
                row.status_to,
                row.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                row.reason,
            ]
            result.append(log_entry)
        return result

    def get_job_logs(self, cromwell_run_ids):
        """
        Get all job logs for a run, with entries organized by job id and ordered by status_from.
        :param cromwell_run_ids: List of Cromwell run IDs (main and any subworkflows)
        :type cromwell_run_ids: list
        :return: job ids and ordered list of log entries (status from, status to, timestamp, reason)
        :rtype: dict
        """
        jobs = {}
        job_logs = self._get_job_logs(cromwell_run_ids)
        for (
            cromwell_run_id,
            cromwell_job_id,
            status_from,
            status_to,
            timestamp,
            reason,
        ) in job_logs:
            cromwell_job_id = str(cromwell_job_id)
            if cromwell_job_id not in jobs:
                jobs[cromwell_job_id] = {}
            index = job_status_value[status_from]
            jobs[cromwell_job_id][index] = [status_from, status_to, timestamp, reason]

        result = {}
        for cromwell_job_id in jobs:
            result[cromwell_job_id] = []
            for index in sorted(jobs[cromwell_job_id].keys()):
                row = jobs[cromwell_job_id][index]
                result[cromwell_job_id].append(row)
        return result

    def _get_cromwell_run_id(self, run_id: int):
        """Get the cromwell_run_id associated with the jaws run_id from the RDb.
        A run may not have a cromwell_run_id if it hasn't been procesed by Cromwell yet.
        An exception is raised if the run is not found in the db.
        :param run_id: JAWS Run ID
        :type run_id: int
        :return: cromwell_run_id UUID
        :rtype: str
        """
        try:
            run = self.session.query(Run).filter_by(id=run_id).one_or_none()
        except Exception as error:
            raise TaskLogDbError(
                f"Task log service was unable to query the db to get the cromwell_run_id for run {run_id}: {error}"
            )
        if run:
            return run.cromwell_run_id  # may be None
        else:
            raise TaskLogRunNotFoundError(f"The run {run_id} was not found")

    def _get_task_summary(self, cromwell_run_id: str):
        """Retrieve all tasks from Cromwell metadata for a run.
        :param cromwell_run_id: Cromwell's UUID for a run
        :type cromwell_run_id: str
        :return: table of tasks
        :rtype: list
        """
        try:
            metadata = self.cromwell.get_metadata(cromwell_run_id)
        except Exception as error:
            err_msg = f"The task log service was unable to retrieve run metadata from Cromwell: {error}"
            self.logger.error(err_msg)
            raise (err_msg)
        summary = metadata.task_summary()
        return summary

    def _get_cromwell_run_ids_from_task_summary(self, summary):
        """Given cromwell metadata task summary table, extract list of cromwell_run_ids"""
        cromwell_run_ids = set()
        for row in summary:
            cromwell_run_id = row[0]
            cromwell_run_ids.add(cromwell_run_id)
        return list(cromwell_run_ids)

    def get_task_log(self, run_id: int):
        """Retrieve complete task log for a run.  This adds task_name and attempt information to the job log."""
        main_cromwell_run_id = self._get_cromwell_run_id(run_id)
        if not main_cromwell_run_id:
            raise TaskLogRunNotFoundError(
                f"Run {run_id} does not have a Cromwell run id (yet)"
            )
        tasks = self._get_task_summary(main_cromwell_run_id)
        all_cromwell_run_ids = self._get_cromwell_run_ids_from_task_summary(tasks)
        all_job_logs = self.get_job_logs(all_cromwell_run_ids)

        tasks_and_logs = []
        for (cromwell_run_id, task_name, attempt, cromwell_job_id) in tasks:
            cromwell_job_id = str(cromwell_job_id)
            if cromwell_job_id not in all_job_logs:
                tasks_and_logs.append(
                    [
                        cromwell_run_id,
                        task_name,
                        attempt,
                        cromwell_job_id,
                        "",
                        "",
                        "",
                        "",
                    ]
                )
                continue
            job_logs = all_job_logs[cromwell_job_id]
            for (status_from, status_to, timestamp, reason) in job_logs:
                tasks_and_logs.append(
                    [
                        cromwell_run_id,
                        task_name,
                        attempt,
                        cromwell_job_id,
                        status_from,
                        status_to,
                        timestamp,
                        reason,
                    ]
                )
        return tasks_and_logs

    def get_task_status(self, run_id: int):
        """
        Retrieve the current status of each task by filtering the log to include only the latest state per task.
        """
        tasks_and_logs = self.get_task_log(run_id)
        tasks_and_last_states = []
        last_job_id = 0
        for row in tasks_and_logs:
            job_id = row[3]
            if job_id == last_job_id:
                tasks_and_last_states[-1] = row
            else:
                tasks_and_last_states.append(row)
                last_job_id = job_id
        return tasks_and_last_states


def get_run_status(session, run_id: int) -> str:
    """
    Retrieve the status of all tasks associated with a run and determine the run's state, with respect to tasks.
    If there are no tasks, return None.
    If any tasks are running/complete, return "running"; "queued" otherwise.

    :param run_id: JAWS run ID
    :type run_id: int
    :return: the status of run, as far as tasks are concerned (None, "queued", or "running")
    :rtype: str
    """
    task_log = TaskLog(session)
    tasks = task_log.get_task_status(run_id)
    if len(tasks) == 0:
        return None
    max_task_status_value = 0
    for task in tasks:
        cromwell_run_id, task_name, attempt, cromwell_job_id, status_from, status_to, timestamp, reason = task
        max_task_status_value = max(max_task_status_value, job_status_value[status_to])
    return "queued" if max_task_status_value < 3 else "running"
