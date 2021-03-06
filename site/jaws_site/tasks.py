"""
This class represents the task-logs which are stored in a rdb and referenced by cromwell_run_id instead of jaws run_id.
"""

import logging
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from dateutil import parser
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
    Extensive use of lazy loading is used.
    """

    def __init__(self, session, **kwargs):
        self.cromwell = cromwell.Cromwell(config.conf.get("CROMWELL", "url"))
        self.session = session
        self._run_id = None
        self._cromwell_run_id = None
        self._job_logs = None
        self._cromwell_task_summary = None
        self._cromwell_job_summary = None
        self._task_status = None
        self._task_log = None
        self._cached_tasks = None
        self._task_summary = None

        if "cromwell_run_id" in kwargs:
            self._cromwell_run_id = kwargs["cromwell_run_id"]
        elif "run_id" in kwargs:
            self._run_id = kwargs["run_id"]
        else:
            raise ValueError("cromwell_run_id or run_id required")

    def cromwell_run_id(self):
        if not self._cromwell_run_id:
            self._get_cromwell_run_id()
        return self._cromwell_run_id

    def _get_cromwell_run_id(self):
        """Get the cromwell_run_id associated with the jaws run_id from the RDb.
        A run may not have a cromwell_run_id if it hasn't been procesed by Cromwell yet.
        An exception is raised if the run is not found in the db.
        :return: cromwell_run_id UUID
        :rtype: str
        """
        assert self._run_id
        try:
            run = self.session.query(Run).filter_by(id=self._run_id).one_or_none()
        except SQLAlchemyError as error:
            raise TaskLogDbError(
                f"Task log service unable to query db for cromwell_run_id of run {self._run_id}: {error}"
            )
        if run:
            self._cromwell_run_id = run.cromwell_run_id  # may be None
            return self._cromwell_run_id
        else:
            raise TaskLogRunNotFoundError(f"The run {self._run_id} was not found")

    def run_id(self):
        if not self._run_id:
            self._get_jaws_run_id()
        return self._run_id

    def _get_jaws_run_id(self):
        """Get the jaws run_id associated with the cromwell_run_id from the RDb.
        An exception is raised if the run is not found in the db.
        """
        assert self._cromwell_run_id
        try:
            run = (
                self.session.query(Run)
                .filter_by(cromwell_run_id=self._cromwell_run_id)
                .one_or_none()
            )
        except SQLAlchemyError as error:
            raise TaskLogDbError(
                f"Task log service was unable to query the db to get the run_id for {self._cromwell_run_id}: {error}"
            )
        if run:
            self._run_id = run.id
        else:
            raise TaskLogRunNotFoundError(
                f"The cromwell run {self._cromwell_run_id} was not found"
            )

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
        cromwell_job_id: int,
        status_from: str,
        status_to: str,
        timestamp: str,
        reason: str = None,
    ):
        """Save a task state transition log"""
        cromwell_run_id = self.cromwell_run_id()
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

    def _select_job_logs(self):
        """Select job logs from db"""
        cromwell_run_id = self.cromwell_run_id()
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
        unsorted_job_logs = []
        for row in table:
            log_entry = [
                row.cromwell_job_id,
                row.status_from,
                row.status_to,
                row.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                row.reason,
            ]
            unsorted_job_logs.append(log_entry)
        return unsorted_job_logs

    def _get_job_logs(self):
        """Get and sort all logs associated with the query cromwell run id."""
        unsorted_job_logs = self._select_job_logs()

        # sort job logs
        jobs = {}
        for (
            cromwell_job_id,
            status_from,
            status_to,
            timestamp,
            reason,
        ) in unsorted_job_logs:
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
        self._job_logs = sorted_logs

    def job_logs(self):
        """
        Get all job state transition logs for a run, with entries organized by job id and ordered by status_from.
        """
        if not self._job_logs:
            self._get_job_logs()
        return self._job_logs

    def cromwell_task_summary(self):
        """Lazy loading of cromwell task summary"""
        if not self._cromwell_task_summary:
            self._get_cromwell_task_summary()
        return self._cromwell_task_summary

    def _get_cromwell_task_summary(self):
        """Retrieve all tasks from Cromwell metadata for a run."""
        cromwell_run_id = self.cromwell_run_id()
        try:
            metadata = self.cromwell.get_metadata(cromwell_run_id)
        except CromwellError as error:
            err_msg = f"The task log service was unable to retrieve run metadata from Cromwell: {error}"
            self.logger.error(err_msg)
            raise (err_msg)
        self._cromwell_task_summary = metadata.task_summary()

    def cromwell_job_summary(self):
        """Return task info, organized by cromwell_job_id."""
        if not self._cromwell_job_summary:
            cromwell_task_summary = self.cromwell_task_summary()
            cromwell_job_summary = {}
            for task_name, cromwell_job_id, cached, max_time in cromwell_task_summary:
                if cromwell_job_id:
                    cromwell_job_id = str(cromwell_job_id)
                    cromwell_job_summary[cromwell_job_id] = [task_name, max_time]
            self._cromwell_job_summary = cromwell_job_summary
        return self._cromwell_job_summary

    def cached_tasks(self):
        """Retrieve all cached tasks from task summary."""
        if self._cached_tasks:
            return self._cached_tasks

        cromwell_task_summary = self.cromwell_task_summary()
        cached_tasks = []
        for task_name, cromwell_job_id, cached, max_time in cromwell_task_summary:
            if cached:
                cached_tasks.append(task_name)
        self._cached_tasks = cached_tasks
        return self._cached_tasks

    def task_log(self):
        """Retrieve complete task log for a run.  This adds task_name to the job log.
        :return: table of task state transitions for the run, including subworkflows
        :rtype: list
        """
        run_id = self.run_id()
        logger.info(f"Run {run_id}: get task-log")
        if self._task_log:
            return self._task_log

        cromwell_run_id = self.cromwell_run_id()
        if not cromwell_run_id:
            # there is no cromwell_run_id if it hasn't been submitted to Cromwell yet
            # (e.g. still in "uploading" state)
            self._task_log = []
            return []

        # Cromwell metadata contains cromwell_job_id and task_name.
        # Note: There can be a delay between job submission and when the job appears in the
        # Cromwell metadata, so some items may be missing.
        cromwell_job_summary = self.cromwell_job_summary()

        # cached tasks don't have job id
        cached_tasks = self.cached_tasks()
        merged_logs = []
        for task_name in cached_tasks:
            merged_logs.append([task_name, None, True, None, None, None, None])

        # The record of job state transitions is stored in a separate db.
        job_logs = self.job_logs()

        # Combine the task names with the logs to produce the final table,
        # ordered by job_id (i.e. order of computation).
        cached = False  # none of these tasks were cached
        for cromwell_job_id in sorted(job_logs.keys()):
            state_transitions = job_logs[cromwell_job_id]
            # default values are required because a job many not appear in the Cromwell metadata immediately
            task_name = "<pending>"
            max_time = None
            if cromwell_job_id in cromwell_job_summary:
                task_name, max_time = cromwell_job_summary[cromwell_job_id]
            for (
                status_from,
                status_to,
                timestamp,
                reason,
            ) in state_transitions:
                merged_logs.append(
                    [
                        task_name,
                        cromwell_job_id,
                        cached,
                        status_from,
                        status_to,
                        timestamp,
                        reason,
                    ]
                )
        self._task_log = merged_logs
        return merged_logs

    def task_summary(self):
        """Retrieve complete task summary for a run.
        :return: foreach task_name, return list of is-cached, queue-time, run-time, result
        :rtype: dict
        """
        run_id = self.run_id()
        logger.info(f"Run {run_id}: get task-summary")
        if self._task_summary:
            return self._task_summary

        task_log = self.task_log()
        task_timestamps = {}
        # task_timestamps value fields are:
        # [0] : cromwell_job_id (string)
        # [1] : cached (bool)
        # [2] : queued (datetime)
        # [3] : running (datetime)
        # [4] : completed (datetime)
        # [5] : result (string)
        for row in task_log:
            (
                task_name,
                cromwell_job_id,
                cached,
                status_from,
                status_to,
                timestamp,
                comment,
            ) = row
            if task_name not in task_timestamps:
                task_timestamps[task_name] = [
                    cromwell_job_id,
                    cached,
                    None,
                    None,
                    None,
                    None,
                ]
            if status_to == "queued":
                task_timestamps[task_name][2] = timestamp
            elif status_to == "running":
                task_timestamps[task_name][3] = timestamp
            elif status_to == "success":
                task_timestamps[task_name][4] = timestamp
                task_timestamps[task_name][5] = "success"
            elif status_to == "failed":
                task_timestamps[task_name][4] = timestamp
                task_timestamps[task_name][5] = "failed"
            if cached is True:
                task_timestamps[task_name][5] = "success"

        cromwell_job_summary = self.cromwell_job_summary()

        task_summary = []
        # task_summary rows' fields are:
        # [0] : task_name (string)
        # [1] : cromwell_job_id (string)
        # [2] : cached (boolean)
        # [3] : result (string)
        # [4] : queued (datetime)
        # [5] : queue-wait (string)
        # [6] : wallclock (string)
        # [7] : max_time (string)
        for task_name, row in task_timestamps.items():
            (cromwell_job_id, cached, queued, running, completed, result) = row
            summary_row = [
                task_name,
                cromwell_job_id,
                cached,
                result,
                queued,
                None,
                None,
                None,
            ]
            if queued and running:
                delta = parser.parse(running) - parser.parse(queued)
                summary_row[5] = str(delta)
            if running and completed:
                delta = parser.parse(completed) - parser.parse(running)
                summary_row[6] = str(delta)
            if cromwell_job_id in cromwell_job_summary:
                summary_row[7] = cromwell_job_summary[cromwell_job_id][1]
            task_summary.append(summary_row)

        self._task_summary = task_summary
        return self._task_summary

    def task_status(self):
        """
        Retrieve the current status of each task by filtering the log to include only the latest state per task.
        """
        run_id = self.run_id()
        logger.info(f"Run {run_id}: Get task-status")
        if self._task_status:
            return self._task_status

        task_log = self.task_log()
        task_status = []
        last_job_id = 0
        for (
            task_name,
            job_id,
            cached,
            status_from,
            status_to,
            timestamp,
            reason,
        ) in task_log:
            # task-status excludes status_from
            row = [task_name, job_id, cached, status_to, timestamp, reason]
            if job_id == last_job_id:
                # state transitions are ordered, so we just keep the last state transition
                task_status[-1] = row
            else:
                task_status.append(row)
                last_job_id = job_id
        self._task_status = task_status
        return self._task_status


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
    task_log = TaskLog(session, run_id=run_id)
    task_status = task_log.task_status()
    if len(task_status) == 0:
        return None
    max_task_status_value = 0
    for task_name, job_id, cached, status, timestamp, reason in task_status:
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
