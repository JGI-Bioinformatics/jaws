from datetime import datetime, timezone
import json
import pika
from time import sleep


DATETIME_FMT = "%Y-%m-%d %H:%M:%S"
SLEEP_STEP = 30


class Sender:
    def __init__(self, config, logger, queue: str):
        self.config = config
        self.logger = logger
        self.queue = queue
        self.sleep_sec = 0

    def send(self, operation: str, params: dict) -> None:
        """
        Login to RabbitMQ server, end a params, and disconnect.
        """
        if type(params) is not dict:
            self.logger.error(f"Discarding invalid params; expected dict: {params}")
            return

        if "timestamp" not in params:
            # timestamps are always in UTC
            params["timestamp"] = datetime.now(timezone.utc).strftime(DATETIME_FMT)

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


class Consumer:
    def __init__(self, config, session, logger, queue: str, operations: dict) -> None:
        """
        Pull messages from the RabbitMQ queue and process them, indefinately.
        """
        self.config = config
        self.session = session
        self.logger = logger
        self.queue = queue
        self.operations = operations
        self.sleep_sec = 0

    def process(self, params: dict) -> bool:
        """
        Process a message by initializing the appropriate object and calling it's update method.
        :param params: The parameters encoded in the message.
        :ptype params: dict
        :return: True if params was processed; False otherwise.
        :rtype: bool
        """
        # datetimes are always stored as UTC and converted to localtime by the query methods
        if "timestamp" not in params:
            params["timestamp"] = datetime.now(timezone.utc)
        elif type(params["timestamp"]) is str:
            params["timestamp"] = datetime.strptime(params["timestamp"], DATETIME_FMT)

        # multiple operations are supported
        operation = params.get("operation", None)
        if operation is None:
            self.logger.error(
                f"Discarding message as missing required 'type': {params}"
            )
        elif operation in self.operations:
            for param in self.operations[operation]["required_params"]:
                if param not in params:
                    self.logger.error(
                        f"Discarding message for {operation} as missing {param}"
                    )
            return self.operations[operation]["function"](self.session, self.logger, params)
        else:
            self.logger.error(f"Discarding message with unknown operation: {operation}")
        return True

    def consume(self):
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
            # close old db session to get a fresh one
            self.session.remove()

            payload = body.decode()
            messages = json.loads(payload)
            if type(messages) is not list:
                self.logger.error(f"Discard message with invalid format: {payload}")
                return True

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
        self.logger.debug("Waiting for messages")
        channel.start_consuming()
