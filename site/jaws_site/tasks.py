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
memory_re = re.compile(r"^\s*(\d*\.?\d*)\s*(\w+)\s*$")


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
                self.session.query(models.Tasks.task_dir)
                .filter(models.Tasks.cromwell_run_id == self.cromwell_run_id)
                .order_by(models.Tasks.id)
            )
        except NoResultFound:
            # it's possible for a newly created run to not have any task-log messages yet
            return None
        except SQLAlchemyError as error:
            self.logger.error(f"Unable to select task_logs: {error}")
            raise TaskDbError(error)
        else:
            return query

    def _select_table(self):
        """
        Select all rows associated with the parent cromwell_run_id; this shall include subworkflows.
        :return: table
        :rtype: list
        """
        try:
            query = (
                self.session.query(models.Tasks)
                .filter(models.Tasks.cromwell_run_id == self.cromwell_run_id)
                .order_by(models.Tasks.id)
            )
        except NoResultFound:
            return []
        except SQLAlchemyError as error:
            self.logger.error(f"Unable to select task_logs: {error}")
            raise TaskDbError(error)

        table = []
        for row in query:
            table.append(
                [
                    row.task_dir,
                    row.status,
                    row.queue_start,
                    row.run_start,
                    row.run_end,
                    row.queue_minutes,
                    row.run_minutes,
                    row.cached,
                    row.name,
                    row.req_cpu,
                    row.req_mem_gb,
                    row.req_minutes,
                ]
            )
        return table

    def _select_all_cpu_minutes(self):
        """
        Select the cpus reserved and the minutes consumed for every task.
        :return: table
        :rtype: list
        """
        try:
            query = (
                self.session.query(models.Tasks)
                .filter(models.Tasks.cromwell_run_id == self.cromwell_run_id)
                .order_by(models.Tasks.id)
            )
        except NoResultFound:
            return []
        except SQLAlchemyError as error:
            self.logger.error(f"Unable to select task_logs: {error}")
            raise TaskDbError(error)
        cpu_minutes = []
        for row in query:
            if row.req_cpu and row.run_minutes:
                cpu_minutes.append([row.req_cpu, row.run_minutes])
        return cpu_minutes

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

    def table(self, **kwargs):
        """
        Convert the timestamps to local timezone strings and return with header.
        """
        rows = self._select_table()
        if "local_tz" in kwargs:
            self.set_local_tz(kwargs.get("local_tz"))
        new_table = []
        for row in rows:
            (
                task_dir,
                status,
                queue_start,
                run_start,
                run_end,
                queue_minutes,
                run_minutes,
                cached,
                name,
                req_cpu,
                req_mem_gb,
                req_minutes,
            ) = row
            queued_str = self._utc_to_local_str(queue_start)
            run_start_str = self._utc_to_local_str(run_start)
            run_end_str = self._utc_to_local_str(run_end)
            new_table.append(
                [
                    task_dir,
                    status,
                    queued_str,
                    run_start_str,
                    run_end_str,
                    queue_minutes,
                    run_minutes,
                    cached,
                    name,
                    req_cpu,
                    req_mem_gb,
                    req_minutes,
                ]
            )
        result = {
            "header": [
                "TASK_DIR",
                "STATUS",
                "QUEUE_START",
                "RUN_START",
                "RUN_END",
                "QUEUE_MIN",
                "RUN_MIN",
                "CACHED",
                "TASK_NAME",
                "REQ_CPU",
                "REQ_GB",
                "REQ_MIN",
            ],
            "data": new_table,
        }
        return result

    def _select_num_started(self):
        try:
            query = (
                self.session.query(models.Tasks)
                .filter(models.Tasks.cromwell_run_id == self.cromwell_run_id)
                .filter(models.Tasks.run_start.isnot(None))
                .all()
            )
        except NoResultFound:
            return 0
        except SQLAlchemyError as error:
            self.logger.error(f"Unable to select task_logs: {error}")
            raise TaskDbError(error)
        else:
            return len(query)

    def did_run_start(self):
        """
        Check if any task has started running by checking the task log.
        """
        return self._select_num_started() > 0

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
            gb = float(m.group(1))
            units = m.group(2).upper()
            if units in ("KB", "KIB"):
                gb = gb / 1024**2
            elif units in ("MB", "MIB"):
                gb = gb / 1024
            elif units in ("GB", "GIB"):
                pass
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

    def _insert_cached_tasks(self, task_summary: dict) -> None:
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
                log_entry = models.Tasks(
                    task_dir=task_dir,
                    name=summary["name"],
                    cromwell_run_id=self.cromwell_run_id,
                    status=status,
                    cached=True,
                    req_cpu=int(summary["req_cpu"]),
                    req_mem_gb=self.memory_gb(summary["req_memory"]),
                    req_minutes=self.time_minutes(summary["req_time"]),
                )
                self.session.add(log_entry)

    def add_metadata(self, summary: dict):
        """
        :param summary: JAWS Cromwell Metadata.task_log_summary() output
        :ptype summary: dict
        """
        data = self._select_rows()
        savepoint = self.session.begin_nested()
        self._insert_cached_tasks(summary)
        if data is not None:
            for rec in data:
                task_dir = rec.task_dir
                if task_dir in summary:
                    info = summary[task_dir]
                    self.logger.debug(f"GOT: {info}")
                    # rec.cached = bool(info["cached"])
                    # rec.name = info["name"]
                    # rec.req_cpu = int(info["requested_cpu"])
                    # rec.req_mem_gb = self.memory_gb(info["requested_memory"])
                    # rec.req_minutes = self.time_minutes(info["requested_time"])
                    # status = info["execution_status"]
                    # if status == "Done":
                    #     rec.status = "succeeded"
                    # elif status == "Failed":
                    #     rec.status = "failed"
                    # elif status == "Aborted":
                    #     rec.status = "cancelled"
        try:
            self.session.commit()
        except SQLAlchemyError as error:
            savepoint.rollback()
            self.logger.exception(
                f"Unable to update Tasks with metadata: {error}"
            )
            raise TaskDbError(f"Unable to update with metadata: {error}")

    def cpu_hours(self) -> float:
        """
        Calculate the total resources consumed (i.e. reserved).
        :return: Total CPU*hours of all tasks.
        :rtype: float
        """
        rows = self._select_all_cpu_minutes()
        cpu_minutes = 0
        for row in rows:
            cpu_minutes = cpu_minutes + row[0] * row[1]
        return round(cpu_minutes / 60, 1)
