"""
This class represents the task-logs which are stored in a rdb and referenced by cromwell_run_id instead of jaws run_id.
"""

import logging
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from jaws_site.models import Job_Log, Run
from jaws_site import config
from jaws_site import cromwell
from jaws_site.cromwell import CromwellError


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
    """This is raised on db errors"""

    pass


class TaskLog:
    """
    The task log is a record of state transitions for tasks.
    """

    def __init__(self, session):
        self.cromwell = cromwell.Cromwell(config.conf.get("CROMWELL", "url"))
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
        except SQLAlchemyError as error:
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
        except SQLAlchemyError as error:
            raise TaskLogError(
                f"Error initializing job_log {cromwell_job_id}:{status_to}: {error}"
            )
        self._save_job_log(job_log)

    def _get_job_logs(self, cromwell_run_id):
        """Get all logs associated with the query cromwell run id.
        :param cromwell_run_id: Cromwell-assigned UUID for the run
        :type cromwell_run_id: str
        :return: table of job logs
        :rtype: list
        """
        try:
            table = (
                self.session.query(Job_Log)
                .filter(Job_Log.cromwell_run_id == cromwell_run_id)
                .all()
            )
        except SQLAlchemyError as error:
            raise TaskLogError(
                f"SQLAlchemy error retrieving task logs from db for run {cromwell_run_id}: {error}"
            )
        except SQLAlchemyError as error:
            raise TaskLogError(
                f"Unknown error retrieving task logs from db for run {cromwell_run_id}: {error}"
            )
        result = []
        for row in table:
            log_entry = [
                row.cromwell_job_id,
                row.status_from,
                row.status_to,
                row.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                row.reason,
            ]
            result.append(log_entry)
        return result

    def get_job_logs(self, cromwell_run_id):
        """
        Get all job state transition logs for a run, with entries organized by job id and ordered by status_from.
        :param cromwell_run_id: Cromwell-assigned UUID for the run
        :type cromwell_run_id: str
        :return: job ids and ordered list of log entries
        :rtype: dict
        """
        jobs = {}
        job_logs = self._get_job_logs(cromwell_run_id)
        for (
            cromwell_job_id,
            status_from,
            status_to,
            timestamp,
            reason,
        ) in job_logs:
            cromwell_job_id = str(cromwell_job_id)
            if cromwell_job_id not in jobs:
                jobs[cromwell_job_id] = {}
            # since job state transition logs come from the backend service, we cannot guarantee every transition is
            # represented in the log, so save in a dict instead of a list, using ordinal value of status_from as key
            index = job_status_value[status_from]
            jobs[cromwell_job_id][index] = [
                status_from,
                status_to,
                timestamp,
                reason,
            ]

        # for each job, create sorted list of it's state transitions
        sorted_logs = {}
        for cromwell_job_id in jobs:
            sorted_logs[cromwell_job_id] = []
            for index in sorted(jobs[cromwell_job_id].keys()):
                row = jobs[cromwell_job_id][index]
                sorted_logs[cromwell_job_id].append(row)
        return sorted_logs

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
        except SQLAlchemyError as error:
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
        except CromwellError as error:
            err_msg = f"The task log service was unable to retrieve run metadata from Cromwell: {error}"
            self.logger.error(err_msg)
            raise (err_msg)
        return metadata.task_summary()

    def get_task_info(self, tasks: list):
        """Retrieve all jobs from Cromwell metadata for a run and reorganize by cromwell_job_id.
        :param tasks: task summary
        :type tasks: list
        :return: cromwell_job_id and task_name
        :rtype: dict
        """
        task_info = {}
        for task_name, cromwell_job_id, cached, run_time in tasks:
            if cromwell_job_id:
                cromwell_job_id = str(cromwell_job_id)
                task_info[cromwell_job_id] = [task_name, run_time]
        return task_info

    def get_cached_tasks(self, tasks: list):
        """Retrieve all cached tasks from task summary.
        :param tasks: task summary
        :type tasks: list
        :return: names to cached tasks
        :rtype: list
        """
        cached_tasks = []
        for task_name, cromwell_job_id, cached, run_time in tasks:
            if cached:
                cached_tasks.append(task_name)
        return cached_tasks

    def get_task_log(self, run_id: int):
        """Retrieve complete task log for a run.  This adds task_name to the job log.
        :param run_id: JAWS Run ID
        :type run_id: int
        :return: table of task state transitions for the run, including subworkflows
        :rtype: list
        """
        logger.info(f"Run {run_id}: get task-log")

        # job logs and cromwell task logs are referenced by cromwell_run_id, so we need to look this up in db
        cromwell_run_id = self._get_cromwell_run_id(run_id)
        if not cromwell_run_id:
            # there is no cromwell_run_id if it hasn't been submitted to Cromwell yet
            # (e.g. still in "uploading" state)
            return []

        # Cromwell metadata contains cromwell_job_id and task_name.
        # Note: There can be a delay between job submission and when the job appears in the
        # Cromwell metadata, so some items may be missing.
        task_summary = self._get_task_summary(cromwell_run_id)

        task_info = self.get_task_info(task_summary)

        # cached tasks don't have job id
        cached_tasks = self.get_cached_tasks(task_summary)
        merged_logs = []
        for task_name in cached_tasks:
            merged_logs.append([task_name, None, None, None, None, "Cached call"])

        # The record of job state transitions is stored in a separate db.
        job_logs = self.get_job_logs(cromwell_run_id)

        # Combine the task names with the logs to produce the final table,
        # ordered by job_id (i.e. order of computation).
        for cromwell_job_id in sorted(job_logs.keys()):
            state_transitions = job_logs[cromwell_job_id]
            # default values are required because a job many not appear in the Cromwell metadata immediately
            task_name = "<pending>"
            run_time = None
            if cromwell_job_id in task_info:
                task_name, run_time = task_info[cromwell_job_id]
            for (
                status_from,
                status_to,
                timestamp,
                reason,
            ) in state_transitions:
                if status_to == "success" and run_time:
                    if reason:
                        reason = f"{reason}; run_time={run_time}"
                    else:
                        reason = f"run_time={run_time}"
                merged_logs.append(
                    [
                        task_name,
                        cromwell_job_id,
                        status_from,
                        status_to,
                        timestamp,
                        reason,
                    ]
                )
        return merged_logs

    def get_task_status(self, run_id: int):
        """
        Retrieve the current status of each task by filtering the log to include only the latest state per task.
        """
        merged_logs = self.get_task_log(run_id)
        tasks_and_last_states = []
        last_job_id = 0
        for task_name, job_id, status_from, status_to, timestamp, reason in merged_logs:
            # task-status excludes status_from
            row = [task_name, job_id, status_to, timestamp, reason]
            if job_id == last_job_id:
                # state transitions are ordered, so we just keep the last state transition
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
    for task_name, job_id, status, timestamp, reason in tasks:
        if status and status in job_status_value:
            max_task_status_value = max(max_task_status_value, job_status_value[status])
        else:
            logger.warn(
                f"Run {run_id} get task status, job {job_id} has unknown status: {status}"
            )
    return "queued" if max_task_status_value < 4 else "running"


def _select_task_log_error_messages(session, cromwell_job_id: str) -> list:
    """
    Get all non-NULL strings from 'reason' column of job logs.
    :param cromwell_job_id: Cromwell's "jobId" field
    :type cromwell_job_ids: str
    :return: list of error messages
    :rtype: list
    """
    try:
        table = (
            session.query(Job_Log)
            .filter_by(cromwell_job_id=cromwell_job_id)
            .filter(Job_Log.reason.isnot(None))
            .filter(Job_Log.reason != "")
            .all()
        )
    except SQLAlchemyError as error:
        raise TaskLogError(
            f"SQLAlchemy error retrieving task logs from db for run {cromwell_job_id}: {error}"
        )
    except Exception as error:
        raise TaskLogError(
            f"Unknown error retrieving task logs from db for run {cromwell_job_id}: {error}"
        )
    result = []
    for row in table:
        result.append(row.reason)
    return result


def get_task_log_error_messages(session, cromwell_job_id):
    """
    Query task-log table any return any error messages found in the optional "reason" column.
    """
    try:
        error_messages = _select_task_log_error_messages(session, cromwell_job_id)
    except Exception as error:
        logger.error(
            f"Error retrieving task log error messages for {cromwell_job_id}: {error}"
        )
        return None
    else:
        return error_messages
