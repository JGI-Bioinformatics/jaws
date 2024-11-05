import logging
from datetime import datetime
from typing import Optional

from jaws_common.exceptions import JawsDbUnavailableError
from sqlalchemy.exc import IntegrityError, OperationalError, SQLAlchemyError
from jaws_site import models


class TaskLogger:
    def __init__(self, config, session, **kwargs) -> None:
        """
        The task logger receives RabbitMQ messages and saves them in the RDb.
        :param config: Configuration parameters
        :ptype config: jaws_site.config
        :param self.session: Database self.session handle
        :ptype self.session: sqlalchemy.orm.self.sessionmaker
        """
        self.config = config
        self.session = session
        self.logger = kwargs.get("logger", None)
        if self.logger is None:
            self.logger = logging.getLogger(__package__)

    def save(self, params: dict) -> bool:
        """
        Save a task log message to the database.
        :param params: Task log message parameters
        :ptype params: dict
        :return: True if the message was saved successfully, False otherwise
        :rtype: bool
        """
        try:
            timestamp = params.get("timestamp")
            params["timestamp"] = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
            status = params.get("status")
            if status not in ["queued", "running", "done"]:
                raise ValueError(f"Invalid status: {status}")

            if status == "queued":
                self.logger.debug(
                    f"Inserting a new task log with queueud status: {str(params)}"
                )
                self._insert(**params)
            else:
                self.logger.debug(f"Updating a tasks's status: {str(params)}")
                self._update(**params)
            return True
        except Exception as e:
            self.logger.error(f"Error saving task log message: {e}")
            return False

    def _insert(self, **kwargs) -> None:
        """
        Insert a new task log entry.

        :param kwargs: Task log entry data
        :type kwargs: dict
        """
        required_keys = ["cromwell_run_id", "task_dir", "status"]
        if not all(key in kwargs for key in required_keys):
            raise ValueError(
                "Missing required keys: {}".format(", ".join(required_keys))
            )

        cromwell_run_id = kwargs.get("cromwell_run_id")
        task_dir = kwargs.get("task_dir")
        timestamp = kwargs.get("timestamp")

        if not self._task_exists(cromwell_run_id, task_dir):
            try:
                log_entry = models.Tasks(
                    cromwell_run_id=cromwell_run_id,
                    task_dir=task_dir,
                    status="queued",
                    queue_start=timestamp,
                )
                self.session.add(log_entry)
                self.session.commit()
            except OperationalError as error:
                # this is the only case in which we would not want to ack the message
                self.logger.error(f"Unable to connect to db: {error}")
                raise JawsDbUnavailableError(str(error))
            except IntegrityError as error:
                self.session.rollback()
                self.logger.error(f"Invalid task-log message, {kwargs}: {error}")
                raise
            except SQLAlchemyError as error:
                self.session.rollback()
                self.logger.exception(
                    f"Failed to insert task log for Task {cromwell_run_id} {task_dir}: {error}"
                )
                raise

    def _task_exists(self, cromwell_run_id: str, task_dir: str) -> bool:
        """
        Check if a task already exists in the database.
        :param cromwell_run_id: Cromwell run ID
        :param task_dir: Task directory
        :return: True if the task exists, False otherwise
        """
        q = self.session.query(models.Tasks).filter(
            models.Tasks.cromwell_run_id == cromwell_run_id,
            models.Tasks.task_dir == task_dir,
        )
        return self.session.query(q.exists()).scalar()

    @staticmethod
    def delta_minutes(start, end) -> Optional[int]:
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

    def _update(self, **kwargs) -> None:
        required_keys = ["cromwell_run_id", "task_dir", "status"]
        if not all(key in kwargs for key in required_keys):
            raise ValueError(
                "Missing required keys: {}".format(", ".join(required_keys))
            )

        cromwell_run_id = kwargs.get("cromwell_run_id")
        task_dir = kwargs.get("task_dir")
        status = kwargs.get("status")
        timestamp = kwargs.get("timestamp")
        return_code = kwargs.get("return_code")
        input_dir_size = kwargs.get("input_dir_size")
        output_dir_size = kwargs.get("output_dir_size")

        row = self._get_task_row(cromwell_run_id, task_dir)

        if row is None:
            self.logger.warning(f"Task {cromwell_run_id} {task_dir} not found!")
        else:
            try:
                if status == "running":
                    row.run_start = timestamp
                    row.status = "running"
                    row.queue_minutes = self.delta_minutes(
                        row.queue_start, row.run_start
                    )
                else:
                    row.run_end = timestamp
                    row.status = "done"
                    row.run_minutes = self.delta_minutes(row.run_start, row.run_end)
                    row.return_code = return_code
                    row.input_dir_size = input_dir_size
                    row.output_dir_size = output_dir_size
                self.session.commit()
            except OperationalError as error:
                # this is the only case in which we would not want to ack the message
                self.session.rollback()
                self.logger.error(f"Unable to connect to db: {error}")
                raise JawsDbUnavailableError(str(error))
            except SQLAlchemyError as error:
                self.session.rollback()
                self.logger.exception(
                    f"Unable to update Tasks {cromwell_run_id} {task_dir}: {error}"
                )
                raise

    def _get_task_row(
        self, cromwell_run_id: str, task_dir: str
    ) -> Optional[models.Tasks]:
        """
        Retrieve the task row for the given cromwell_run_id and task_dir.
        :param cromwell_run_id: Cromwell run ID
        :param task_dir: Task directory
        :return: The task row if found, None otherwise
        """
        try:
            return (
                self.session.query(models.Tasks)
                .filter(
                    models.Tasks.cromwell_run_id == cromwell_run_id,
                    models.Tasks.task_dir == task_dir,
                )
                .order_by(models.Tasks.id.desc())
                .first()
            )
        except OperationalError as error:
            self.logger.error(f"Unable to connect to db: {error}")
            raise JawsDbUnavailableError(str(error))
        except Exception as error:
            self.logger.error(
                f"Unexpected error; unable to select task log {cromwell_run_id} {task_dir}: {error}"
            )
            raise
