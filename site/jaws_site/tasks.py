import logging
from datetime import datetime, timezone
import pytz
import re
from sqlalchemy import update
from sqlalchemy.exc import IntegrityError, SQLAlchemyError, OperationalError
from sqlalchemy.orm.exc import NoResultFound
from jaws_site import models


DEFAULT_TZ = "America/Los_Angeles"
DATETIME_FMT = "%Y-%m-%d %H:%M:%S"
time_re = re.compile(r"^\s*(\d+):(\d+):(\d+)\s*$")
memory_re = re.compile(r"^\s*(\d*\.?\d*)\s*(\w+)\s*$")
DEFAULT_CPU = 1
DEFAULT_MEM_GB = 5


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
        self.data = None

    def set_local_tz(self, local_tz=DEFAULT_TZ):
        self.local_tz = local_tz
        self.local_tz_obj = pytz.timezone(local_tz)

    def select(self):
        """
        Select all rows associated with the parent cromwell_run_id; this shall include subworkflows.
        :return: table
        :rtype: list
        """
        if self.data is not None:
            return self.data
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
                    row.id,
                    row.task_dir,
                    row.job_id,
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
        self.data = table
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
        if "local_tz" in kwargs:
            self.set_local_tz(kwargs.get("local_tz"))
        new_table = []
        rows = self.select()
        for row in rows:
            (
                row_id,
                task_dir,
                job_id,
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
            cpu_hours = None
            if req_cpu is not None and run_minutes is not None:
                cpu_hours = round(req_cpu * run_minutes / 60, 3)
            new_table.append(
                [
                    task_dir,
                    job_id,
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
                    cpu_hours,
                ]
            )
        result = {
            "header": [
                "TASK_DIR",
                "JOB_ID",
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
                "CPU_HRS",
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
    def memory_gb(memory_str: str) -> float:
        """
        Convert memory to gigabytes.
        :param memory_str: number and units (e.g. "5 GB", "512 MB", "0.5 TB")
        :ptype memory_str: str
        :return: memory in gigabytes
        :rtype: float
        """
        if memory_str is not None:
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
        return None

    @staticmethod
    def int_or_none(x: any) -> int:
        """
        :param x: a variable, presumably an integer
        :ptype x: any
        :return: integer value or None
        :rtype: int
        """
        if x is None:
            return None
        elif type(x) is int:
            return x
        elif type(x) is str and x.isdigit():
            return int(x)
        elif type(x) is float:
            return int(round(x, 0))
        else:
            return None

    def _insert_cached_tasks(self, summary: dict) -> None:
        """
        Cached tasks don't appear in the call-log, so we'll copy them from the metadata for completeness.
        They shall not include cpu-hours data so won't affect resource calculations.
        New rows are inserted but not committed as this is part of a larger transaction; see: add_metadata().
        :param summary: list of task metadata, produced by JAWS Cromwell class
        :ptype summary: list
        """
        for task_dir, info in summary.items():
            if info["cached"] is True:
                status = info["execution_status"]
                if status == "Done":
                    status = "succeeded"
                elif status == "Failed":
                    status = "failed"
                elif status == "Aborted":
                    status = "cancelled"
                req_cpu = self.int_or_none(info.get("requested_cpu", DEFAULT_CPU))
                req_mem_gb = self.memory_gb(
                    info.get("requested_memory", DEFAULT_MEM_GB)
                )
                req_minutes = self.int_or_none(
                    info.get("requested_runime_minutes", None)
                )
                log_entry = models.Tasks(
                    task_dir=task_dir,
                    name=info["name"],
                    cromwell_run_id=self.cromwell_run_id,
                    status=status,
                    cached=True,
                    req_cpu=req_cpu,
                    req_mem_gb=req_mem_gb,
                    req_minutes=req_minutes,
                )
                self.session.add(log_entry)

    def prepare_metadata(self, summary: dict) -> list:
        """
        Prepare bulk update list
        """
        updates = []
        rows = self.select()
        for row in rows:
            (
                row_id,
                task_dir,
                job_id,
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
            if task_dir not in summary:
                continue
            status = "cancelled"
            if summary[task_dir]["execution_status"] == "Done":
                status = "succeeded"
            elif summary[task_dir]["execution_status"] == "Failed":
                status = "failed"
            req_cpu = self.int_or_none(
                summary[task_dir].get("requested_cpu", DEFAULT_CPU)
            )
            req_mem_gb = self.memory_gb(
                summary[task_dir].get("requested_memory", DEFAULT_MEM_GB)
            )
            req_minutes = self.int_or_none(
                summary[task_dir].get("requested_runtime_minutes", None)
            )
            update = {
                "id": row_id,
                "job_id": job_id,
                "status": status,
                "cached": bool(summary[task_dir]["cached"]),
                "name": summary[task_dir]["name"],
                "req_cpu": req_cpu,
                "req_mem_gb": req_mem_gb,
                "req_minutes": req_minutes,
            }
            updates.append(update)
        return updates

    def add_metadata(self, summary: dict):
        """
        :param summary: JAWS Cromwell Metadata.task_log_summary() output
        :ptype summary: dict
        """
        savepoint = self.session.begin_nested()
        self._insert_cached_tasks(summary)
        updates = self.prepare_metadata(summary)
        self.logger.debug(f"UPDATES: {updates}")
        try:
            self.session.commit()
            self.session.execute(update(models.Tasks), updates)
        except SQLAlchemyError as error:
            savepoint.rollback()
            self.logger.exception(f"Unable to update Tasks with metadata: {error}")
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
        return round(cpu_minutes / 60, 2)


def save_task_log(session, logger, **kwargs) -> bool:
    required_params = ["cromwell_run_id", "task_dir", "status", "timestamp"]
    for param in required_params:
        if param not in kwargs:
            logger.error(f"Invalid task log, missing 'param': {kwargs}")
            return True
    if kwargs["status"] == "queued":
        return insert_task_log(session, logger, **kwargs)
    else:
        return update_task_log(session, logger, **kwargs)


def insert_task_log(session, logger, **kwargs) -> bool:
    """
    Insert new row.
    :return: True if saved or discarded; False otherwise.
    :rtype: bool
    """
    cromwell_run_id = kwargs.get("cromwell_run_id")
    task_dir = kwargs.get("task_dir")
    timestamp = kwargs.get("timestamp")

    try:
        log_entry = models.Tasks(
            cromwell_run_id=cromwell_run_id,
            task_dir=task_dir,
            status="queued",
            queue_start=timestamp,
        )
        session.add(log_entry)
        session.commit()
    except OperationalError as error:
        # this is the only case in which we would not want to ack the message
        logger.error(f"Unable to connect to db: {error}")
        return False
    except IntegrityError as error:
        session.rollback()
        logger.error(f"Invalid task-log message, {kwargs}: {error}")
    except SQLAlchemyError as error:
        session.rollback()
        logger.exception(
            f"Failed to insert task log for Task {cromwell_run_id} {task_dir}: {error}"
        )
    return True


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


def update_task_log(session, logger, **kwargs) -> bool:
    """
    Update a task's log to record it's new state.
    :return: True if saved or discarded; False otherwise.
    :rtype: bool
    """
    cromwell_run_id = kwargs.get("cromwell_run_id")
    task_dir = kwargs.get("task_dir")
    status = kwargs.get("status")
    timestamp = kwargs.get("timestamp")

    # select row for this task
    try:
        row = (
            session.query(models.Tasks)
            .filter(
                models.Tasks.cromwell_run_id == cromwell_run_id,
                models.Tasks.task_dir == task_dir,
            )
            .one_or_none()
        )
    except OperationalError as error:
        # this is the only case in which we would not want to ack the message
        logger.error(f"Unable to connect to db: {error}")
        return False
    except Exception as error:
        err_msg = f"Unexpected error; unable to select task log {cromwell_run_id} {task_dir}: {error}"
        logger.error(err_msg)
        return True

    if row is None:
        # this can only happen if the original "queued" message was lost
        # TODO: should we insert a record with the current time instead?
        logger.error(f"Task {cromwell_run_id} {task_dir} not found!")
        return True

    try:
        if status == "running":
            row.run_start = timestamp
            row.status = "running"
            row.queue_minutes = delta_minutes(row.queue_start, row.run_start)
        else:
            row.run_end = timestamp
            row.status = "done"
            row.run_minutes = delta_minutes(row.run_start, row.run_end)
        session.commit()
    except OperationalError as error:
        # this is the only case in which we would not want to ack the message
        session.rollback()
        logger.error(f"Unable to connect to db: {error}")
        return False
    except SQLAlchemyError as error:
        session.rollback()
        logger.exception(
            f"Unable to update Tasks {cromwell_run_id} {task_dir}: {error}"
        )
    return True
