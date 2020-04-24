"""
A Scalable and threaded Consumer that will automatically re-connect on failure.
"""
import logging
import threading
import time
import json
import urllib.parse
import amqpstorm

from jaws_site import config
from jaws_site.dispatch import dispatch


class Consumer(object):

    def __init__(self, rpc_queue):
        """Initialize Consumer object

        :param rpc_queue: The name of the queue from which to retrieve messages.
        :type rpc_queue: str
        :return:
        """
        self.logger = logging.getLogger(__package__)
        self.rpc_queue = rpc_queue
        self.channel = None
        self.active = False

    def start(self, connection):
        """Start the consumer"""
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
        """Stop the consumer"""
        if self.channel:
            self.channel.close()

    def __call__(self, message):
        """Process the RPC Payload.
        :param message: A JSON-RPC2 encoded request
        :type message: str
        :return:
        """
        # DECODE JSON-RPC2 REQUEST
        request = json.loads(message.body)
        method = request["method"]
        params = request["params"]

        # GET RESPONSE FROM DISPATCHER
        response_dict = dispatch(method, params)

        # VALIDATE RESPONSE
        response_dict["jsonrpc"] = "2.0"
        if "error" in response_dict:
            if type(response_dict["error"]) is not dict:
                self.logger.error(f'Invalid response_dict: {response_dict}')
            elif "code" not in response_dict["error"]:
                self.logger.error("Invalid response_dict; missing error code")
            elif "message" not in response_dict["error"]:
                self.logger.error("Invalid response_dict; missing error message")
            elif "result" in response_dict:
                self.logger.error("Invalid response_dict; result not allowed if error")
        elif "result" not in response_dict:
            self.logger.error(f'Invalid response_dict: {response_dict}')

        # INIT RESPONSE MESSAGE OBJECT
        properties = {"correlation_id": message.correlation_id}
        response = amqpstorm.Message.create(
            message.channel, json.dumps(response_dict), properties
        )

        # DELIVER RESPONSE
        response.publish(message.reply_to)
        message.ack()


class RpcServer(object):

    def __init__(self) -> None:
        self.logger = logging.getLogger(__package__)
        self.hostname = config.conf.get("AMQP", "host")
        self.vhost = config.conf.get("AMQP", "vhost")
        self.user = config.conf.get("AMQP", "user")
        self.password = config.conf.get("AMQP", "password")
        self.rpc_queue = config.conf.get("AMQP", "queue")
        self.max_retries = int(config.conf.get("RPC", "max_retries"))
        number_of_consumers = int(config.conf.get("RPC", "num_threads"))
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
                self.logger.warning(str(e))
                self.stop_consumers(len(self.consumers))
                self.create_connection()
            time.sleep(1)

    def increase_consumers_by(self, num) -> None:
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
                    uri = f'amqp://{self.user}:' \
                          f'{urllib.parse.quote_plus(self.password)}@' \
                          f'{self.hostname}:5672/{self.vhost}?heartbeat=6'
                    return amqpstorm.UriConnection(uri)
                else:
                    return amqpstorm.Connection(self.hostname, self.user,
                                                self.password)
            except amqpstorm.AMQPConnectionError as e:
                self.logger.warning(str(e))
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
        for _ in num:
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
