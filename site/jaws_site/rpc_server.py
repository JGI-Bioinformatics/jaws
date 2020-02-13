"""
A Scalable and threaded Consumer that will automatically re-connect on failure.
"""
import logging
import threading
import time
import json
import os
import urllib.parse
import amqpstorm
from amqpstorm import Connection
from amqpstorm import Message

from jaws_site import dispatch

logging.basicConfig(level=logging.INFO)

LOGGER = logging.getLogger()

DEBUG = True if "JAWS_DEBUG" in os.environ else False

class RpcServer(object):

    def __init__(self, config):
        self.config = config
        self.hostname = config["AMQP"]["host"]
        self.vhost = config["AMQP"]["vhost"]
        self.username = config["AMQP"]["user"]
        self.password = config["AMQP"]["password"]
        self.rpc_queue = config["AMQP"]["queue"]
        self.number_of_consumers = int(config["RPC"]["num_threads"])
        self.max_retries = int(config["RPC"]["max_retries"])
        self._connection = None
        self._consumers = []
        self._stopped = threading.Event()

    def start_server(self):
        """Start the RPC Server.
        :return:
        """
        self._stopped.clear()
        if not self._connection or self._connection.is_closed:
            self._create_connection()
        while not self._stopped.is_set():
            try:
                # Check our connection for errors.
                self._connection.check_for_errors()
                self._update_consumers()
            except amqpstorm.AMQPError as why:
                # If an error occurs, re-connect and let update_consumers
                # re-open the channels.
                LOGGER.warning(why)
                self._stop_consumers()
                self._create_connection()
            time.sleep(1)

    def increase_consumers(self):
        """Add one more consumer.
        :return:
        """
        if self.number_of_consumers <= 20:
            self.number_of_consumers += 1

    def decrease_consumers(self):
        """Stop one consumer.
        :return:
        """
        if self.number_of_consumers > 0:
            self.number_of_consumers -= 1

    def stop(self):
        """Stop all consumers.
        :return:
        """
        while self._consumers:
            consumer = self._consumers.pop()
            consumer.stop()
        self._stopped.set()
        self._connection.close()

    def _create_connection(self):
        """Create a connection.
        :return:
        """
        attempts = 0
        while True:
            attempts += 1
            if self._stopped.is_set():
                break
            try:
                #self._connection = Connection(self.hostname, self.username, self.password)
                if self.vhost:
                    uri = "amqp://" + self.username + ":" + urllib.parse.quote_plus(self.password) + "@" + self.hostname + ":5672/" + self.vhost + "?heartbeat=60"
                    self._connection = amqpstorm.UriConnection(uri)
                else:
                    self._connection = amqpstorm.Connection(self.hostname, self.username, self.password)
                break
            except amqpstorm.AMQPError as why:
                LOGGER.warning(why)
                if self.max_retries and attempts > self.max_retries:
                    raise Exception('max number of retries reached')
                time.sleep(min(attempts * 2, 30))
            except KeyboardInterrupt:
                break

    def _update_consumers(self):
        """Update Consumers.
            - Add more if requested.
            - Make sure the consumers are healthy.
            - Remove excess consumers.
        :return:
        """
        # Do we need to start more consumers.
        consumer_to_start = \
            min(max(self.number_of_consumers - len(self._consumers), 0), 2)
        for _ in range(consumer_to_start):
            consumer = Consumer(self.rpc_queue, self.config)
            self._start_consumer(consumer)
            self._consumers.append(consumer)

        # Check that all our consumers are active.
        for consumer in self._consumers:
            if consumer.active:
                continue
            self._start_consumer(consumer)
            break

        # Do we have any overflow of consumers.
        self._stop_consumers(self.number_of_consumers)

    def _stop_consumers(self, number_of_consumers=0):
        """Stop a specific number of consumers.
        :param number_of_consumers:
        :return:
        """
        while len(self._consumers) > number_of_consumers:
            consumer = self._consumers.pop()
            consumer.stop()

    def _start_consumer(self, consumer):
        """Start a consumer as a new Thread.
        :param Consumer consumer:
        :return:
        """
        thread = threading.Thread(target=consumer.start, args=(self._connection,))
        thread.daemon = True
        thread.start()


class Consumer(object):
    def __init__(self, rpc_queue, config):
        self.rpc_queue = rpc_queue
        self.channel = None
        self.active = False
        self.dispatcher = dispatch.Dispatcher()

    def start(self, connection):
        if DEBUG: print("Starting consumer")
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

    def __call__(self, message):
        """Process the RPC Payload.
        :param Message message:
        :return:
        """
        corr_id = message.correlation_id

        # DECODE JSON-RPC2 REQUEST
        request = json.loads(message.body)
        method = request["method"]
        params = request["params"]

        if DEBUG: print("Dispatching request for method %s, corr_id %s" % (method, corr_id))

        # GET RESPONSE FROM DISPATCHER
        response_dict = self.dispatcher.dispatch(method, params)

        # VALIDATE RESPONSE
        response_dict["jsonrpc"] = "2.0"
        if "error" in response_dict:
            if type(response_dict["error"]) is not dict: sys.exit("Invalid response_dict: %s" % (response_dict,))
            if "code" not in response_dict["error"]: sys.exit("Invalid response_dict; missing error code")
            if "message" not in response_dict["error"]: sys.exit("Invalid response_dict; missing error message")
            if "result" in response_dict: sys.exit("Invalid response_dict; result not allowed if error")
        elif "result" not in response_dict: sys.exit("Invalid response_dict: %s" % (response_dict,))

        # INIT RESPONSE MESSAGE OBJECT
        properties = { 'correlation_id': message.correlation_id }
        response = Message.create(message.channel, json.dumps(response_dict), properties)

        # DELIVER RESPONSE
        if DEBUG: print("Publishing response for %s" % (corr_id,))
        response.publish(message.reply_to)
        message.ack()
