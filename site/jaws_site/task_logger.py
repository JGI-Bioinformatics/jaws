import json
import logging
import pika
from datetime import datetime
from sqlalchemy.exc import IntegrityError, SQLAlchemyError, OperationalError
from jaws_site import models


DATETIME_FMT = "%Y-%m-%d %H:%M:%S"


class TaskLoggerDbError(Exception):
    pass


class TaskLogger:
    def __init__(self, config, session, logger=None) -> None:
        """
        The task logger receives RabbitMQ messages and saves them in the RDb.
        :param config: Configuration parameters
        :ptype config: jaws_site.config
        :param self.session: Database self.session handle
        :ptype self.session: sqlalchemy.orm.self.sessionmaker
        """
        self.config = config
        self.session = session
        if logger is None:
            self.logger = logging.getLogger(__package__)
        else:
            self.logger = logger

    def save(self, message: str) -> bool:
        params = json.loads(message)
        timestamp = params.get("timestamp")
        params["timestamp"] = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
        status = params.get("status")

        if status == "queued":
            return self._insert(**params)
        elif status in ("running", "done"):
            return self._update(**params)
        else:
            self.logger.error(f"Invalid log status: {status}")
            return True

    def _insert(self, **kwargs) -> bool:
        """
        Insert new row.
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
        except OperationalError as error:
            # this is the only case in which we would not want to ack the message
            self.logger.error(f"Unable to connect to db: {error}")
            return False
        except IntegrityError as error:
            self.session.rollback()
            self.logger.error(f"Invalid task-log message, {kwargs}: {error}")
        except SQLAlchemyError as error:
            self.session.rollback()
            self.logger.exception(
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

    def _update(self, **kwargs) -> bool:
        cromwell_run_id = kwargs.get("cromwell_run_id")
        task_dir = kwargs.get("task_dir")
        status = kwargs.get("status")
        timestamp = kwargs.get("timestamp")

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
            return False
        except Exception as error:
            err_msg = f"Unexpected error; unable to select task log {cromwell_run_id} {task_dir}: {error}"
            self.logger.error(err_msg)
            return True

        if row is None:
            # this can only happen if the queued message was lost
            self.logger.error(f"Task {cromwell_run_id} {task_dir} not found!")
            return True

        try:
            if status == "running":
                row.run_start = timestamp
                row.status = "running"
                row.queue_minutes = self.delta_minutes(row.queue_start, row.run_start)
            else:
                row.run_end = timestamp
                row.status = "done"
                row.run_minutes = self.delta_minutes(row.run_start, row.run_end)
            self.session.commit()
        except OperationalError as error:
            # this is the only case in which we would not want to ack the message
            self.session.rollback()
            self.logger.error(f"Unable to connect to db: {error}")
            return False
        except SQLAlchemyError as error:
            self.session.rollback()
            self.logger.exception(
                f"Unable to update Tasks {cromwell_run_id} {task_dir}: {error}"
            )
        return True

    def receive_messages(self):
        """
        Receive and consume task-log messages indefinately.
        """
        host = self.config.get("RMQ", "host")
        port = self.config.get("RMQ", "port")
        vhost = self.config.get("RMQ", "vhost")
        user = self.config.get("RMQ", "user")
        password = self.config.get("RMQ", "password")
        ttl = int(
            self.config.get("RMQ", "ttl", 3600000)
        )  # must match value in sender's self.config
        site_id = self.config.get("SITE", "id")
        queue = f"{site_id}_tasks"

        def callback(ch, method, properties, body):
            message = body.decode()
            self.logger.debug(str(message))
            if self.save(message) is True:
                ch.basic_ack(delivery_tag=method.delivery_tag)

        credentials = pika.PlainCredentials(user, password)
        connection = pika.BlockingConnection(
            pika.ConnectionParameters(host, port, vhost, credentials)
        )
        channel = connection.channel()

        try:
            channel.queue_declare(
                queue=queue, durable=True, arguments={"x-message-ttl": ttl}
            )
        except Exception as error:
            self.logger.warning(f"Unable to declare RMQ queue: {error}")
            # channel.queue_delete(queue=queue)
            # connection.close()
            raise

        channel.basic_consume(queue=queue, on_message_callback=callback, auto_ack=False)
        self.logger.debug("Waiting for task-log messages")
        channel.start_consuming()