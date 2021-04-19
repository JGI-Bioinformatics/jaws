"""
This class represents the task-logs which are stored in a rdb and referenced by cromwell_run_id instead of jaws run_id.
"""

import logging
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from jaws_site.models import Job_Log, Run
from jaws_site import config
from jaws_site import cromwell


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
        Get all job state transition logs for a run, with entries organized by job id and ordered by status_from.
        :param cromwell_run_ids: List of Cromwell run IDs (main and any subworkflows)
        :type cromwell_run_ids: list
        :return: job ids and ordered list of log entries
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
            # since job state transition logs come from the backend service, we cannot guarantee every transition is
            # represented in the log, so save in a dict instead of a list, using ordinal value of status_from as key
            index = job_status_value[status_from]
            jobs[cromwell_job_id][index] = [
                cromwell_run_id,
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
        except Exception as error:
            raise TaskLogDbError(
                f"Task log service was unable to query the db to get the cromwell_run_id for run {run_id}: {error}"
            )
        if run:
            return run.cromwell_run_id  # may be None
        else:
            raise TaskLogRunNotFoundError(f"The run {run_id} was not found")

    def _get_cromwell_task_summary(self, cromwell_run_id: str):
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
        return metadata.task_summary()

    def get_job_metadata(self, cromwell_run_id: str):
        """Retrieve all jobs from Cromwell metadata for a run and reorganize by cromwell_job_id.
        :param cromwell_run_id: Cromwell's UUID for the run
        :type cromwell_run_id: str
        :return: cromwell_job_id and task metadata
        :rtype: dict
        """
        tasks = self._get_cromwell_task_summary(cromwell_run_id)
        jobs = {}
        for (cromwell_run_id, task_name, attempt, cromwell_job_id) in tasks:
            cromwell_job_id = str(cromwell_job_id)
            jobs[cromwell_job_id] = [cromwell_run_id, task_name, attempt]
        return jobs

    def _get_cromwell_run_ids_from_job_metadata(self, jobs):
        """Given jobs, return list (without duplicates) of cromwell_run_ids.
        :param jobs: cromwell_job_id and associated metadata
        :type jobs: dict
        :return: unique list of cromwell_run_ids
        :rtype: list
        """
        cromwell_run_ids = set()
        for cromwell_job_id in jobs.keys():
            cromwell_run_id = jobs[cromwell_job_id][0]
            cromwell_run_ids.add(cromwell_run_id)
        return list(cromwell_run_ids)

    def get_task_log(self, run_id: int):
        """Retrieve complete task log for a run.  This adds task_name and attempt information to the job log.
        :param run_id: JAWS Run ID
        :type run_id: int
        :return: table of task state transitions for the run, including subworkflows
        :rtype: list
        """
        # job logs and cromwell task logs are referenced by cromwell_run_id, so we need to look this up in db
        main_cromwell_run_id = self._get_cromwell_run_id(run_id)
        if not main_cromwell_run_id:
            # run exists but hasn't been submitted to Cromwell yet (e.g. uploading)
            return []

        # Cromwell metadata contains cromwell_job_id, task_name, and attempt, but does not
        # contain detailed state transitions.  There can be a delay between job submission
        # and when the job appears in the Cromwell metadata, so some items may be missing.
        job_metadata = self.get_job_metadata(main_cromwell_run_id)

        # While we have the cromwell_run_id of the main workflow in the Runs db,
        # if a run has subworkflows, we need to get their cromwell_run_ids from the
        # Cromwell metadata.
        cromwell_run_ids = self._get_cromwell_run_ids_from_job_metadata(job_metadata)

        # Select all logs for the main run and any subworkflows from the Task_Log table.
        # The db contains all job state transitions and is current, but the records do not
        # have "task_name" or "attempt" metadata (only Cromwell metadata has those fields)
        job_logs = self.get_job_logs(cromwell_run_ids)

        # Combine the job metadata (which include task names) with the logs (i.e. state
        # transitions) to produce the final, complete record, which is returned in a table.
        merged_logs = []
        for cromwell_job_id in sorted(job_logs.keys()):
            state_transitions = job_logs[cromwell_job_id]
            # default values are required because a job many not appear in the Cromwell
            # metadata immediately
            task_name = "<pending>"
            attempt = "?"
            if cromwell_job_id in job_metadata:
                cromwell_run_id, task_name, attempt = job_metadata[cromwell_job_id]
            for (
                cromwell_run_id,
                status_from,
                status_to,
                timestamp,
                reason,
            ) in state_transitions:
                merged_logs.append(
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
        return merged_logs

    def get_task_status(self, run_id: int):
        """
        Retrieve the current status of each task by filtering the log to include only the latest state per task.
        """
        merged_logs = self.get_task_log(run_id)
        tasks_and_last_states = []
        last_job_id = 0
        for row in merged_logs:
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
        (
            cromwell_run_id,
            task_name,
            attempt,
            cromwell_job_id,
            status_from,
            status_to,
            timestamp,
            reason,
        ) = task
        max_task_status_value = max(max_task_status_value, job_status_value[status_to])
    return "queued" if max_task_status_value < 3 else "running"


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
            .filter(Job_Log.reason.isnot(""))
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
    return _select_task_log_error_messages(session, cromwell_job_id)
