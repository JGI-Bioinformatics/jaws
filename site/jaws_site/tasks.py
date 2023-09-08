import logging
from datetime import datetime, timezone
import pytz
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from jaws_site import models


DEFAULT_TZ = "America/Los_Angeles"
DATETIME_FMT = "%Y-%m-%d %H:%M:%S"


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
        Select rows, copy into a table, and calculate queue and run durations.  All timestamps are datetime objects in UTC, durations are timedelta objects.
        """
        table = []
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

        for row in query:
            queue_start = row.queue_start
            run_start = row.run_start
            run_end = row.run_end
            queue_dur = None
            run_dur = None
            if queue_start and run_start:
                queue_dur = run_start - queue_start
            if run_start and run_end:
                run_dur = run_end - run_start

            table.append(
                [
                    row.task_dir,
                    row.status,
                    queue_start,
                    run_start,
                    run_end,
                    row.rc,
                    queue_dur,
                    run_dur,
                ]
            )
        return table

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
        table = []
        for row in self.data:
            (
                task_dir,
                status,
                queue_start,
                run_start,
                run_end,
                rc,
                queue_dur,
                run_dur,
            ) = row
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
                    str(queue_dur),
                    str(run_dur),
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
                "QUEUE_DUR",
                "RUN_DUR",
            ],
            "data": table,
        }
        return result

    def did_run_start(self):
        """
        Check if any task has started running by checking the task log.
        """
        for row in self.data:
            run_start = row[3]
            if run_start is not None:
                return True
        return False
