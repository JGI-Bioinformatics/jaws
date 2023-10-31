#!/usr/bin/env python

from datetime import datetime, timezone
from jaws_site import tasks
import json
import logging
import pika
from time import sleep


DATETIME_FMT = "%Y-%m-%d %H:%M:%S"
SLEEP_STEP = 30


class MessageSender:
    def __init__(self, config, logger=None, **kwargs):
        self.config = config
        self.sleep_sec = 0
        if logger is None:
            self.logger = logging.getLogger(__package__)
        else:
            self.logger = logger
        if "queue" in kwargs:
            self.queue = kwargs.get("queue")
        else:
            site_id = self.config.get("SITE", "id")
            deployment = self.config.get("SITE", "deployment")
            self.queue = f"jaws_{site_id}_{deployment}_logs"

    def send_message(self, message: dict) -> None:
        """
        Login to RabbitMQ server, end a message, and disconnect.
        :param message: the message to send
        :ptype message: dict
        """
        if type(message) is not dict:
            self.logger.error(f"Discarding invalid message; expected dict: {message}")
            return

        if "timestamp" not in message:
            # timestamps are always in UTC
            message["timestamp"] = datetime.now(timezone.utc).strftime(DATETIME_FMT)

        host = self.config.get("RMQ", "host")
        port = self.config.get("RMQ", "port")
        vhost = self.config.get("RMQ", "vhost")
        user = self.config.get("RMQ", "user")
        password = self.config.get("RMQ", "password")
        # both queue and ttl in both sender and receiver must be identical
        ttl = int(self.config.get("RMQ", "ttl"))

        credentials = pika.PlainCredentials(user, password)
        connection = pika.BlockingConnection(
            pika.ConnectionParameters(host, port, vhost, credentials)
        )
        channel = connection.channel()
        arguments = {"x-message-ttl": ttl}
        channel.queue_declare(queue=self.queue, durable=True, arguments=arguments)
        channel.basic_publish(
            exchange="",
            routing_key=self.queue,
            body=json.dumps(message),
            properties=pika.BasicProperties(
                delivery_mode=pika.spec.PERSISTENT_DELIVERY_MODE
            ),
        )
        connection.close()


class MessageReceiver:
    def __init__(self, config, session, logger=None, **kwargs) -> None:
        """
        Pull messages from the RabbitMQ queue and process them, indefinately.
        :param config: Configuration parameters
        :ptype config: jaws_site.config
        :param self.session: Database self.session handle
        :ptype self.session: sqlalchemy.orm.self.sessionmaker
        """
        self.config = config
        self.session = session
        self.sleep_sec = 0
        if logger is None:
            self.logger = logging.getLogger(__package__)
        else:
            self.logger = logger
        if "queue" in kwargs:
            self.queue = kwargs.get("queue")
        else:
            site_id = self.config.get("SITE", "id")
            deployment = self.config.get("SITE", "deployment")
            self.queue = f"jaws_{site_id}_{deployment}_logs"

    def process(self, params: dict) -> bool:
        """
        Process a message by initializing the appropriate object and calling it's update method.
        :param params: The parameters encoded in the message.
        :ptype params: dict
        :return: True if params was processed; False otherwise.
        :rtype: bool
        """
        # datetimes are always stored as UTC and converted to localtime by the query methods
        if "timestamp" in params:
            # if timestamp (str) was provided, convert to datetime object
            params["timestamp"] = datetime.strptime(params["timestamp"], DATETIME_FMT)
        else:
            params["timestamp"] = datetime.now(timezone.utc)

        # multiple operations are supported
        operation = params.get("operation", None)
        if operation is None:
            self.logger.error(
                f"Discarding message as missing required 'type': {params}"
            )
        elif operation in operations:
            for param in operations[operation]["required_params"]:
                if param not in params:
                    self.logger.error(
                        f"Discarding message for {operation} as missing {param}"
                    )
            return operations[operation]["function"](self.session, self.logger, params)
        else:
            self.logger.error(f"Discarding message with unknown operation: {operation}")
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
        # both queue and ttl in both sender and receiver must be identical
        ttl = int(self.config.get("RMQ", "ttl", 3600000))

        def callback(ch, method, properties, body):
            """
            Decode and process messages.  Acknowledge unless RDb error.
            """
            payload = body.decode()
            messages = json.loads(payload)
            okay = True
            for message in messages:
                if self.process(message) is False:
                    okay = False
                    break
            if okay is True:
                ch.basic_ack(delivery_tag=method.delivery_tag)
                # reset sleep amount after each success
                self.sleep_sec = 0
            else:
                self.sleep_sec = self.sleep_sec + SLEEP_STEP
                self.logger.warn(f"RDb error; sleep for {self.sleep_src}")
                sleep(self.sleep_src)

        credentials = pika.PlainCredentials(user, password)
        connection = pika.BlockingConnection(
            pika.ConnectionParameters(host, port, vhost, credentials)
        )
        channel = connection.channel()

        try:
            channel.queue_declare(
                queue=self.queue, durable=True, arguments={"x-message-ttl": ttl}
            )
        except Exception as error:
            self.logger.warning(f"Unable to declare RMQ queue: {error}")
            # channel.queue_delete(queue=self.queue)
            # connection.close()
            raise

        channel.basic_consume(
            queue=self.queue, on_message_callback=callback, auto_ack=False
        )
        self.logger.debug("Waiting for task-log messages")
        channel.start_consuming()


operations = {
    "task_log": {
        "required_params": ["cromwell_run_id", "task_dir", "status", "timestamp"],
        "function": tasks.save_task_log,
    },
}
