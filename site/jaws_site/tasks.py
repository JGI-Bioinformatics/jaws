import json
import logging
import pika
from datetime import datetime
from dateutil import parser
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from jaws_site import models


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
            table.append(
                [
                    row.execution_dir,
                    row.status,
                    row.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                ]
            )
        self.data = table

    def table(self):
        """
        Calculate queue/run durations and reformat into a table with columns:
        - cromwell_run_id
        - status
        - queued timestamp
        - running timestamp
        - completed timestamp
        - queue duration
        - running duration
        """
        execution_dirs = {}
        for (execution_dir, status, timestamp) in self.data:
            if execution_dir not in execution_dirs:
                execution_dirs[execution_dir] = ["", "", "", ""]
            if status == "queued":
                execution_dirs[execution_dir][1] = timestamp
            elif status == "running":
                execution_dirs[execution_dir][2] = timestamp
            else:
                execution_dirs[execution_dir][3] = timestamp
                execution_dirs[execution_dir][0] = status
        table = []
        for execution_dir in sorted(execution_dirs.keys()):
            # calculate queue and run durations
            row = execution_dirs[execution_dir]
            if row[1] is not None and row[2] is not None:
                delta = parser.parse(row[2]) - parser.parse(row[1])
                row.append(str(delta))
            else:
                row.append(None)
            if row[2] is not None and row[3] is not None:
                delta = parser.parse(row[3]) - parser.parse(row[2])
                row.append(str(delta))
            else:
                row.append(None)
            # add execution dir
            row = (execution_dir, *execution_dirs[execution_dir])
            table.append(row)
        return table

    def did_run_start(self):
        """
        Check if any task has started running by checking the task log.
        """
        rows = self._select_rows()
        for row in rows:
            (cromwell_run_id, execution_dir, status, timestamp) = row
            if status != "queued":
                return True
        return False


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
    ttl = int(config.get("RMQ", "ttl", 3600000))  # must match value in sender's config
    site_id = config.get("SITE", "id")
    queue = f"{site_id}_tasks"

    def _insert_task_log(message: str) -> None:
        params = json.loads(message)
        timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
        cromwell_run_id = params.get("cromwell_run_id", None)
        execution_dir = params.get("execution_dir", None)
        try:
            log_entry = models.Task_Log(
                cromwell_run_id=cromwell_run_id,
                execution_dir=execution_dir,
                status=params.get("status"),
                timestamp=timestamp,
            )
            session.add(log_entry)
            session.commit()
        except IntegrityError as error:
            session.rollback()
            logger.error(f"Invalid task-log message, {message}: {error}")
        except SQLAlchemyError as error:
            session.rollback()
            logger.exception(
                f"Failed to insert task log for Task {execution_dir} ({status}): {error}"
            )

    def callback(ch, method, properties, body):
        message = body.decode()
        logger.debug(f"Received: {message}")
        try:
            _insert_task_log(message)
        except Exception as error:
            logger.error(f"Unable to save task log: {error}")
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
