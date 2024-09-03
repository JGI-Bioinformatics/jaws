import logging
from datetime import datetime

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
        Save or update task information based on the provided parameters.

        :param params: Dictionary containing task information
        :return: True if the operation was successful, False otherwise
        :raises ValueError: If the status is invalid
        """
        try:
            timestamp = params.get("timestamp")
            params["timestamp"] = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
            
            status = params.get("status")
            if status not in ["queued", "running", "done"]:
                raise ValueError(f"Invalid status: {status}")

            if status == "queued":
                return self._insert(**params)
            else:
                return self._update(**params)
        except ValueError as e:
            self.logger.error(f"Invalid input: {e}")
            raise
        except Exception as e:
            self.logger.exception(f"Unexpected error in save method: {e}")
            return False

    def _insert(self, **kwargs) -> bool:
        """
        Insert new row into the Tasks table.
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
            self.session.add(log_entry)
            self.session.commit()
            self.logger.info(f"Successfully inserted task log for {cromwell_run_id} {task_dir}")
            return True
        except OperationalError as error:
            self.logger.error(f"Database connection error: {error}")
            raise JawsDbUnavailableError(f"Unable to connect to database: {error}") from error
        except IntegrityError as error:
            self.session.rollback()
            self.logger.error(f"Integrity error for task-log message {kwargs}: {error}")
            return False
        except SQLAlchemyError as error:
            self.session.rollback()
            self.logger.exception(f"Failed to insert task log for {cromwell_run_id} {task_dir}: {error}")
            return False
        except Exception as error:
            self.session.rollback()
            self.logger.exception(f"Unexpected error inserting task log for {cromwell_run_id} {task_dir}: {error}")
            return False

    @staticmethod
    def delta_minutes(start: datetime, end: datetime) -> int | None:
        """
        Calculate the difference between two timestamps, rounded to the nearest minute.

        :param start: Start time
        :param end: End time
        :return: Difference in minutes (rounded) or None if either input is None
        """
        if start is None or end is None:
            return None
        
        duration = end - start
        return round(duration.total_seconds() / 60)

    def _update(self, **kwargs) -> bool:
        cromwell_run_id = kwargs.get("cromwell_run_id")
        task_dir = kwargs.get("task_dir")
        status = kwargs.get("status")
        timestamp = kwargs.get("timestamp")
        return_code = kwargs.get("return_code")

        # select row for this task
        try:
            row = (
                self.session.query(models.Tasks)
                .filter(
                    models.Tasks.cromwell_run_id == cromwell_run_id,
                    models.Tasks.task_dir == task_dir,
                )
                .one_or_none()
            )
        except OperationalError as error:
            # this is the only case in which we would not want to ack the message
            self.logger.error(f"Unable to connect to db: {error}")
            raise JawsDbUnavailableError(str(error))
        except Exception as error:
            err_msg = f"Unexpected error; unable to select task log {cromwell_run_id} {task_dir}: {error}"
            self.logger.error(err_msg)
            raise

        if row is None:
            # this can only happen if the queued message was lost
            self.logger.error(f"Task {cromwell_run_id} {task_dir} not found!")
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
