import logging
from datetime import datetime, timezone
import pytz
import re
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from jaws_site import models


DEFAULT_TZ = "America/Los_Angeles"
DATETIME_FMT = "%Y-%m-%d %H:%M:%S"
time_re = re.compile(r"^\s*(\d+):(\d+):(\d+)\s*$")
memory_re = re.compile(r"^\s*(\d+)\s*(\w+)\s*$")


class TaskDbError(Exception):
    pass


class TaskLog:
    """
    Whenever a Run's task changes state (e.g. is submitted to cluster, starts running,
    finishes), a message is sent via RabbitMQ.  This class receives those messages and
    saves them in a persistent RDb and provides them upon request.
    """

    def __init__(self, session, cromwell_run_id, logger=None, **kwargs) -> None:
        self.session = session
        self.cromwell_run_id = cromwell_run_id
        if logger is None:
            self.logger = logging.getLogger(__package__)
        else:
            self.logger = logger
        self.set_local_tz(kwargs.get("local_tz", DEFAULT_TZ))
        self.data = self._select_rows()

    def set_local_tz(self, local_tz=DEFAULT_TZ):
        self.local_tz = local_tz
        self.local_tz_obj = pytz.timezone(local_tz)

    def _select_rows(self):
        """
        Select all rows associated with the parent cromwell_run_id; this shall include subworkflows.
        :return: sqlalchemy query object
        :rtype: obj
        """
        try:
            query = (
                self.session.query(models.Task_Log)
                .filter(models.Task_Log.cromwell_run_id == self.cromwell_run_id)
                .order_by(models.Task_Log.id)
            )
        except NoResultFound:
            # it's possible for a newly created run to not have any task-log messages yet
            return []
        except SQLAlchemyError as error:
            self.logger.error(f"Unable to select task_logs: {error}")
            raise TaskDbError(error)
        else:
            return query

    def _utc_to_local_str(self, timestamp) -> str:
        """Convert UTC time to the local time zone. This should handle daylight savings.
        :param timestamp: a datetime object, UTC timezone
        :ptype timestamp: datetime.datetime
        :return: a formatted string in the local timezone
        :rtype: str
        """
        if timestamp is None:
            return None
        if type(timestamp) is str:
            timestamp = datetime.strptime("2023-04-24 08:00:00", DATETIME_FMT)
        local_timestamp = timestamp.replace(tzinfo=timezone.utc).astimezone(
            tz=self.local_tz_obj
        )
        return local_timestamp.strftime(DATETIME_FMT)

    @staticmethod
    def delta_minutes(start, end) -> int:
        """
        Return the difference between two timestamps, rounded to the nearest minute.
        :param start: start time
        :ptype: datetime.datetime
        :param end: end time
        :ptype end: datetime.datetime
        :return: difference in minutes (rounded)
        :rtype: int
        """
        if start and end:
            duration = end - start
            return round(duration.total_seconds() / 60, 0)
        else:
            return None

    def table(self, **kwargs):
        """
        Convert the timestamps to local timezone strings and return with header.
        """
        if "local_tz" in kwargs:
            self.set_local_tz(kwargs.get("local_tz"))
        table = []
        for row in self.data:
            task_dir = row.task_dir
            status = row.status
            rc = row.rc
            queue_start = row.queue_start
            run_start = row.run_start
            run_end = row.run_end
            queue_minutes = self.delta_minutes(queue_start, run_start)
            run_minutes = self.delta_minutes(run_start, run_end)
            queued_str = self._utc_to_local_str(queue_start)
            run_start_str = self._utc_to_local_str(run_start)
            run_end_str = self._utc_to_local_str(run_end)
            table.append(
                [
                    task_dir,
                    status,
                    queued_str,
                    run_start_str,
                    run_end_str,
                    rc,
                    queue_minutes,
                    run_minutes,
                ]
            )
        result = {
            "header": [
                "TASK",
                "STATUS",
                "QUEUE_START",
                "RUN_START",
                "RUN_END",
                "RC",
                "QUEUE_MINUTES",
                "RUN_MINUTES",
            ],
            "data": table,
        }
        return result

    def did_run_start(self):
        """
        Check if any task has started running by checking the task log.
        """
        for row in self.data:
            if row.run_start is not None:
                return True
        return False

    @staticmethod
    def time_minutes(duration: str) -> int:
        """
        Convert time duration to minutes.
        :param duration: time duration in hh:mm:ss format
        :ptype duration: str
        :return: time in minutes
        :rtype: int
        """
        m = time_re.match(duration)
        if m:
            minutes = (
                float(m.group(1)) * 60 + float(m.group(2)) + float(m.group(3)) / 60
            )
            return round(minutes, 0)
        else:
            return None

    @staticmethod
    def memory_gb(memory_str: str) -> float:
        """
        Convert memory to gigabytes.
        :param memory_str: number and units (e.g. "5 GB", "512 MB", "0.5 TB")
        :ptype memory_str: str
        :return: memory in gigabytes
        :rtype: float
        """
        m = memory_re.match(memory_str)
        if m:
            gb = m.group(1)
            units = m.group(2).upper()
            if units in ("KB", "KIB"):
                gb = gb / 1024**2
            elif units in ("MB", "MIB"):
                gb = gb / 1024
            elif units in ("GB", "GIB"):
                gb = float(gb)
            elif units in ("TB", "TIB"):
                gb = gb * 1024
            else:
                gb = gb / 1024**3
            if gb > 32767:
                # max value for signed SMALLINT
                gb = 32767
            return round(gb, 0)
        else:
            return None

    def _insert_cached_taks(self, task_summary: list) -> None:
        """
        Cached tasks don't appear in the call-log, so we'll copy them from the metadata for completeness.
        They shall not include cpu-hours data so won't affect resource calculations.
        New rows are inserted but not committed as this is part of a larger transaction; see: add_metadata().
        :param task_summary: list of task metadata, produced by JAWS Cromwell class
        :ptype task_summary: list
        """
        for task_dir, summary in task_summary.items():
            if summary["cached"] is True:
                status = summary["execution_status"]
                if status == "Done":
                    status = "succeeded"
                elif status == "Failed":
                    status = "failed"
                elif status == "Aborted":
                    status = "cancelled"
                log_entry = models.Task_Log(
                    task_dir=task_dir,
                    name=summary["name"],
                    cromwell_run_id=self.cromwell_run_id,
                    cromwell_job_id=summary["job_id"],
                    status=status,
                    queue_start=None,
                    cached=True,
                    requested_minutes=self.time_minutes(summary["requested_time"]),
                    requested_cpu=int(summary["requested_cpu"]),
                    requested_gb=self.memory_gb(summary["requested_memory"]),
                )
                self.session.add(log_entry)

    def _update_with_cromwell_metadata(self, task_summary: list) -> None:
        """
        Add fields from the Cromwell metadata:
        - requested cpu
        - requested memory (in GB)
        - requested time (in minutes)
        And update the status using the Cromwell "executionStatus",  which shall differ from return code when a
        task's expected outputs cannot be reaped (e.g. missing file).
        """
        for row in self.data:
            task_dir = row.task_dir
            if task_dir in task_summary:
                summary = task_summary[task_dir]
                row.requested_minutes = self.time_minutes(summary["requested_time"])
                row.requested_cpu = int(summary["requested_cpu"])
                row.requested_gb = self.memory_gb(summary["requested_memory"])
                status = summary["execution_status"]
                if status == "Done":
                    row.status = "succeeded"
                elif status == "Failed":
                    row.status = "failed"
                elif status == "Aborted":
                    row.status = "cancelled"

    def add_metadata(self, metadata):
        """ """
        summary = metadata.task_summary_dict()
        savepoint = self.session.begin_nested()
        try:
            self._insert_cached_tasks(summary)
            self._update_with_cromwell_metadata(summary)
            self.session.commit()
        except SQLAlchemyError as error:
            savepoint.rollback()
            self.logger.exception(
                f"Unable to update Run {self.data.id} Tasks with metadata: {error}"
            )
            raise TaskDbError(f"Unable to update with metadata: {error}")
