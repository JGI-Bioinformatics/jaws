import logging
import re
from datetime import datetime, timezone

import pytz
from sqlalchemy import update
from sqlalchemy.exc import SQLAlchemyError
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
                    row.return_code,
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
                return_code,
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
                    return_code,
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
                "RETURN_CODE",
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
                return_code = info.get("return_code", None)
                log_entry = models.Tasks(
                    task_dir=task_dir,
                    name=info["name"],
                    cromwell_run_id=self.cromwell_run_id,
                    status=status,
                    cached=True,
                    req_cpu=req_cpu,
                    req_mem_gb=req_mem_gb,
                    req_minutes=req_minutes,
                    return_code=return_code,
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
                return_code,
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
            job_id = summary[task_dir].get("job_id", None)

            # metadata's return_code is None mostly as expandSubWorkflows=0
            # return_code_meta = summary[task_dir].get("return_code", None)
            # tasks table's return_code is updated whenever a task is completed
            # from the rc file

            update = {
                "id": row_id,
                "job_id": job_id,
                "status": status,
                "cached": bool(summary[task_dir]["cached"]),
                "name": summary[task_dir]["name"],
                "req_cpu": req_cpu,
                "req_mem_gb": req_mem_gb,
                "req_minutes": req_minutes,
                "return_code": return_code,
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
        self.logger.debug(f"Metadata to update tasks table: {updates}")
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
