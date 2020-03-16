"""
A Scalable and threaded Consumer that will automatically re-connect on failure.
"""
import logging
from typing import List
import threading
import time
import json
import urllib.parse
import amqpstorm
from jaws_site.dispatch import dispatch
from jaws_site.config import JawsConfig

logger = logging.getLogger(__package__)


class Consumer(object):

    rpc_queue: str
    channel: amqpstorm.Channel
    active: bool

    def __init__(self, rpc_queue: str):
        self.rpc_queue = rpc_queue
        self.channel = None
        self.active = False

    def start(self, connection: amqpstorm.Connection):
        self.channel = None
        try:
            self.active = True
            self.channel = connection.channel(rpc_timeout=10)
            self.channel.basic.qos(1)
            self.channel.queue.declare(self.rpc_queue)
            self.channel.basic.consume(self, self.rpc_queue, no_ack=False)
            self.channel.start_consuming()
            if not self.channel.consumer_tags:
                # Only close the channel if there is nothing consuming.
                # This is to allow messages that are still being processed
                # in __call__ to finish processing.
                self.channel.close()
        except amqpstorm.AMQPError:
            pass
        finally:
            self.active = False

    def stop(self):
        if self.channel:
            self.channel.close()

    def __call__(self, message: amqpstorm.Message):
        """Process the RPC Payload.
        :param Message message:
        :return:
        """
        corr_id = message.correlation_id

        # DECODE JSON-RPC2 REQUEST
        request = json.loads(message.body)
        method = request["method"]
        params = request["params"]

        # GET RESPONSE FROM DISPATCHER
        response_dict = dispatch(method, params)

        # VALIDATE RESPONSE
        response_dict["jsonrpc"] = "2.0"
        if "error" in response_dict:
            if not isinstance(response_dict["error"], dict):
                raise InvalidJsonRpcResponse(f'{response_dict["error"]} is not of type dict')
            if "code" not in response_dict["error"]:
                raise InvalidJsonRpcResponse('missing error code')
            if "message" not in response_dict["error"]:
                raise InvalidJsonRpcResponse("missing error message")
            if "result" in response_dict:
                raise InvalidJsonRpcResponse("result not allowed in error")
        elif "result" not in response_dict:
            raise InvalidJsonRpcResponse('missing result in response')

        # INIT RESPONSE MESSAGE OBJECT
        properties = {'correlation_id': corr_id}
        response = amqpstorm.Message.create(message.channel, json.dumps(response_dict), properties)

        # DELIVER RESPONSE
        response.publish(message.reply_to)
        message.ack()


class RpcServer(object):

    hostname: str
    vhost: str
    user: str
    password: str
    rpc_queue: str
    number_of_consumers: int
    max_retries: int
    connection: amqpstorm.Connection
    consumers: List[Consumer]
    stopped: threading.Event

    def __init__(self, conf: JawsConfig) -> None:
        self.hostname = conf.get("AMQP", "host")
        self.vhost = conf.get("AMQP", "vhost")
        self.user = conf.get("AMQP", "user")
        self.password = conf.get("AMQP", "password")
        self.rpc_queue = conf.get("AMQP", "queue")
        self.max_retries = int(conf.get("RPC", "max_retries"))
        number_of_consumers = int(conf.get("RPC", "num_threads"))
        self.consumers = [Consumer(self.rpc_queue) for _ in range(number_of_consumers)]
        self.stopped = threading.Event()
        self.connection = self.create_connection()

    def start_server(self) -> None:
        """Start the RPC Server.
        :return:
        """
        self.stopped.clear()
        while not self.stopped.is_set():
            try:
                # Check our connection for errors.
                self.connection.check_for_errors()
                self.update_consumers()
            except amqpstorm.AMQPConnectionError as e:
                # If an error occurs, re-connect and let update_consumers
                # re-open the channels.
                logger.warning(str(e))
                self.stop_consumers(len(self.consumers))
                self.create_connection()
            time.sleep(1)

    def increase_consumers_by(self, num: int) -> None:
        """Add one more consumer.
        :return: index where we activate consumers
        """
        for _ in range(num):
            consumer = Consumer(self.rpc_queue)
            self.start_consumer(consumer)
            self.consumers.append(consumer)

    def decrease_consumers_by(self, num: int) -> None:
        """Stop one consumer.
        :return:
        """
        if len(self.consumers) - num < 0:
            num = len(self.consumers)
        for _ in range(num):
            consumer = self.consumers.pop()
            consumer.stop()

    def stop(self) -> None:
        """Stop all consumers.
        :return:
        """
        for consumer in self.consumers:
            consumer.stop()
        del self.consumers[:]
        self.stopped.set()
        self.connection.close()

    def create_connection(self) -> amqpstorm.Connection:
        """
        Returns a Connection to a RabbitMQ service

        :return: amqpstorm.Connection
        """
        for _ in range(self.max_retries):
            if self.stopped.is_set():
                break
            try:
                if self.vhost:
                    uri = f'amqp://{self.user}:{urllib.parse.quote_plus(self.password)}@{self.hostname}:5672/{self.vhost}?heartbeat=60' # noqa
                    return amqpstorm.UriConnection(uri)
                else:
                    return amqpstorm.Connection(self.hostname, self.user, self.password)
            except amqpstorm.AMQPConnectionError as e:
                logger.warning(str(e))
                time.sleep(1)
        raise ExceededRetries("Reached the maximum level of retries")

    def update_consumers(self, add_consumer: int = 0):
        """Update Consumers.
            - Add more if requested.
            - Make sure the consumers are healthy.
            - Remove excess consumers.
        :return:
        """
        # Do we need to start more consumers.
        self.increase_consumers_by(add_consumer)

        # Check that all our consumers are active.
        for consumer in self.consumers:
            if not consumer.active:
                self.start_consumer(consumer)

    def stop_consumers(self, num: int):
        """Stop a specific number of consumers.
        :param num number_of_consumers:
        :return:
        """
        for _ in len(num):
            consumer = self.consumers.pop()
            consumer.stop()

    def start_consumer(self, consumer: Consumer):
        """Start a consumer as a new Thread.
        :param Consumer consumer:
        :return:
        """
        thread = threading.Thread(target=consumer.start, args=(self.connection,))
        thread.daemon = True
        thread.start()


class ExceededRetries(Exception):
    pass


class InvalidJsonRpcResponse(Exception):
    pass
