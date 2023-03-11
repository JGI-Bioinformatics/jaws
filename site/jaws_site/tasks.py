import io
import json
import logging
import os
import pika
from datetime import datetime
from sqlalchemy import select
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from jaws_site import models, config


class TaskLog:
    def __init__(self, session, cromwell_run_id, logger=None) -> None:
        self.session = session
        self.cromwell_run_id = cromwell_run_id
        if logger is None:
            self.logger = logging.getLogger(__package__)
        else:
            self.logger = logger
        try:
            self.data = self._select_rows()
        except TaskDbError as error:
            raise error

    def _select_rows(self):
        try:
            stmt = (
                select(
                    models.Task_Log.execution_dir,
                    models.Task_Log.status,
                    models.Task_Log.timestamp,
                )
                .where(
                    models.Task_Log.any(
                        models.Task_Log.cromwell_run_id == self.cromwell_run_id
                    )
                )
                .order_by(models.Task_Log.id)
            )
            rows = self.session.execute(stmt).all()
        except NoResultFound:
            # it's possible for a newly created run to not have any task-log messages yet
            return []
        except SQLAlchemyError as error:
            self.logger.error(f"Unable to select task_logs: {error}")
            raise TaskDbError(error)
        return rows

    def table(self):
        """
        Reformat task log into a table.
        """
        execution_dirs = {}
        for execution_dir, status, timestamp in self.data:
            if not execution_dir in execution_dirs:
                execution_dirs[execution_dir] = ["?", "?", "?", "?"]
                # the "?" should be replaced unless a log message was lost or
                # the task is unfinished
            if status == "queued":
                execution_dirs[execution_dir][0] = timestamp
            elif status == "running":
                execution_dirs[execution_dir][1] = timestamp
            else:
                execution_dirs[execution_dir][2] = timestamp
                execution_dirs[execution_dir][3] = status
        table = []
        for execution_dir in sort(execution_dirs):
            row = (execution_dir, *execution_dirs[execution_dir])
            table.append(row)
        return table


def receive_messages(config, session):
    """
    Receive and consume task-log messages indefinately.
    :param config: Configuration parameters
    :ptype config: jaws_site.config
    :param session: Database session handle
    :ptype session: sqlalchemy.orm.sessionmaker
    """

    logger = logging.getLogger(__package__)

    host = config.get("RMQ", "host")
    port = config.get("RMQ", "port")
    vhost = config.get("RMQ", "vhost")
    user = config.get("RMQ", "user")
    password = config.get("RMQ", "password")
    ttl = int(config.get("RMQ", "ttl"))
    site_id = config.get("SITE", "id")
    queue = f"{site_id}_tasks"

    def _insert_task_log(message: str) -> None:
        (execution_dir, status, timestamp) = json.loads(message)
        cromwell_run_id = "/".split(execution_dir)[0]
        timestamp = datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%S.%f%z")
        try:
            log_entry = models.Task_Log(
                cromwell_run_id=cromwell_run_id,
                execution_dir=execution_dir,
                status=status,
                timestamp=timestamp,
            )
            session.add(log_entry)
            session.commit()
        except SQLAlchemyError as error:
            session.rollback()
            logger.exception(
                f"Failed to insert task log for Task {execution-dir} ({status}): {error}"
            )

    def callback(ch, method, properties, body):
        message = body.decode()
        logger.debug(f"Received: {message}")
        try:
            _insert_task_log(message)
        except:
            logger.error("Unable to save task log")
            raise
        else:
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
        logger.warning(f"Unable to declare RMQ queue: {error}")
        # channel.queue_delete(queue=queue)
        # connection.close()
        raise

    channel.basic_consume(queue=queue, on_message_callback=callback, auto_ack=False)
    logger.debug("Waiting for task-log messages")
    channel.start_consuming()
