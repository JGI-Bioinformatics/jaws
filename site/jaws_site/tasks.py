import logging
from datetime import datetime, timezone
import pytz
from dateutil import parser
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from jaws_site import models


DATETIME_FMT = "%Y-%m-%d %H:%M:%S"


class TaskDbError(Exception):
    pass


class TaskLog:
    """
    Whenever a Run's task changes state (e.g. is submitted to cluster, starts running,
    finishes), a message is sent via RabbitMQ.  This class receives those messages and
    saves them in a persistent RDb and provides them upon request.
    """

    def __init__(self, session, cromwell_run_id, logger=None) -> None:
        self.session = session
        self.cromwell_run_id = cromwell_run_id
        if logger is None:
            self.logger = logging.getLogger(__package__)
        else:
            self.logger = logger
        self.data = self._select_rows()

    def _select_rows(self):
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
            queued = row.queued
            if queued is not None:
                queued = queued.strftime(DATETIME_FMT)
            running = row.running
            if running is not None:
                running = running.strftime(DATETIME_FMT)
            completed = row.completed
            if completed is not None:
                completed = completed.strftime(DATETIME_FMT)
            cancelled = row.cancelled
            if cancelled is not None:
                cancelled = cancelled.strftime(DATETIME_FMT)

            queue_dur = None
            run_dur = None
            if running is not None:
                delta = parser.parse(running) - parser.parse(queued)
                queue_dur = str(delta)
                if completed is not None:
                    delta = parser.parse(completed) - parser.parse(running)
                    run_dur = str(delta)

            table.append(
                [
                    row.task_dir,
                    queued,
                    running,
                    completed,
                    row.rc,
                    cancelled,
                    queue_dur,
                    run_dur,
                ]
            )
        return table

    @staticmethod
    def _utc_to_local(utc_datetime: str, local_tz: str) -> str:
        """Convert UTC time to the local time zone. This should handle daylight savings.
        :param utc_datetime: a string of date and time "2021-07-06 11:15:17".
        :ptype utc_datetime: str
        :param local_tz: name of the local timezone; use server timezone if not specified
        :ptype local_tz: str
        :return: similarly formatted string in the specified local tz
        :rtype: str
        """
        # The timezone can be overwritten with a environmental variable.
        # JAWS_TZ should be set to a timezone in a similar format to 'US/Pacific'
        local_tz_obj = ""
        if local_tz is None:
            local_tz_obj = datetime.now().astimezone().tzinfo
        else:
            local_tz_obj = pytz.timezone(local_tz)
        datetime_obj = datetime.strptime(utc_datetime, DATETIME_FMT)
        local_datetime_obj = datetime_obj.replace(tzinfo=timezone.utc).astimezone(
            tz=local_tz_obj
        )
        return local_datetime_obj.strftime(DATETIME_FMT)

    def table_local_tz(self, local_tz: str) -> list:
        table = []
        for row in self.data:
            (
                task_dir,
                queued,
                running,
                completed,
                rc,
                cancelled,
                queue_dir,
                run_dir,
            ) = row
            queued = self._utc_to_local(queued, local_tz)
            if running is not None:
                running = self._utc_to_local(running, local_tz)
            if completed is not None:
                completed = self._utc_to_local(completed, local_tz)
            if cancelled is not None:
                cancelled = self._utc_to_local(cancelled, local_tz)
            table.append(
                [
                    task_dir,
                    queued,
                    running,
                    completed,
                    rc,
                    cancelled,
                    queue_dir,
                    run_dir,
                ]
            )
        return table

    def table(self, **kwargs):
        """
        Update the times to local, if specified, and return with header.
        """
        local_tz = kwargs.get("local_tz", None)
        table = self.data
        if local_tz is not None:
            table = self.table_local_tz(local_tz)
        result = {
            "header": [
                "TASK",
                "QUEUED",
                "RUNNING",
                "COMPLETED",
                "RC",
                "CANCELLED",
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
            running = row[2]
            if running is not None:
                return True
        return False
