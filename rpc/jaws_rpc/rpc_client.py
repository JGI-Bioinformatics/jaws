"""RPC_Client sends JSON-RPC2 messages to an RPC_Server."""

import threading
import amqpstorm
from amqpstorm import Message
import json
from time import sleep
from jaws_rpc import jsonrpc_utils


DEFAULT_PORT = 5672
DEFAULT_WAIT_INTERVAL = 0.25
DEFAULT_MAX_WAIT = 10
RPC_QUEUE = "rpc"
MESSAGE_TTL = "10"  # expires in seconds


class RPC_Client(object):
    """Asynchronous remote procedure call (RPC) client class."""

    def __init__(self, params):
        """Constructor

        :param params: A dictionary containing configuration parameters.
        :type params: dict
        """
        self.params = params
        self.queue = {}
        self.channel = None
        self.connection = None
        self.callback_queue = None
        self.rpc_queue = RPC_QUEUE
        self.open()
        self.wait_interval = DEFAULT_WAIT_INTERVAL
        if "port" not in self.params:
            self.params["port"] = DEFAULT_PORT
        if "rpc_wait_interval" in self.params:
            self.wait_interval = float(self.params["rpc_wait_interval"])
        else:
            self.wait_interval = DEFAULT_WAIT_INTERVAL
        if "rpc_max_wait" in self.params:
            self.max_wait = float(self.params["rpc_max_wait"])
        else:
            self.max_wait = DEFAULT_MAX_WAIT
        for required_param in ["host", "user", "password"]:
            if required_param not in self.params:
                raise ConfigurationError(f"{required_param} required")

    def open(self):
        """Open connection to RabbitMQ"""
        try:
            self.connection = amqpstorm.Connection(
                self.params["host"],
                self.params["user"],
                self.params["password"],
                int(self.params["port"]),
                virtual_host=self.params["vhost"]
            )
        except Exception as error:
            raise ConnectionError(error)
        self.channel = self.connection.channel()
        self.channel.queue.declare(RPC_QUEUE)
        result = self.channel.queue.declare(exclusive=True)
        self.callback_queue = result["queue"]
        self.channel.basic.consume(
            self._on_response, no_ack=True, queue=self.callback_queue
        )
        self._create_process_thread()

    def _create_process_thread(self):
        """Create a thread responsible for consuming messages in response
         to RPC requests.
        """
        thread = threading.Thread(target=self._process_data_events)
        thread.setDaemon(True)
        thread.start()

    def _process_data_events(self):
        """Process Data Events using the Process Thread."""
        self.channel.start_consuming()

    def _on_response(self, message):
        """On Response store the message with the correlation id in a local
         dictionary.
        """
        self.queue[message.correlation_id] = message.body

    def send_request(self, method, params={}):
        """Format and send a JSON-RPC2 request but do not wait for a response.

        :param method: The RPC method to call; this is not validated by the rpc client.
        :type method: str
        :param params: Any associated parameters to accompany the request (varies by method).
        :type params: dict
        :return: the ID of the RPC request (correlation ID)
        :rtype: int
        """
        # Construct JSON-RPC2 payload string
        request = {"jsonrpc": "2.0", "method": method, "params": params}
        try:
            jsonrpc_utils.is_valid_request(request)
        except Exception:
            raise

        payload = json.dumps(request, default=str)

        # Create the Message object.
        message = Message.create(self.channel, payload, properties={"expiration": MESSAGE_TTL})
        message.reply_to = self.callback_queue

        # Create an entry in our local dictionary, using the automatically
        # generated correlation_id as our key.
        self.queue[message.correlation_id] = None

        # Publish the RPC request.
        message.publish(routing_key=self.rpc_queue)

        # Return the Unique ID used to identify the request.
        return message.correlation_id

    def get_response(self, corr_id):
        """Return the JSON-RPC2 response if ready, None otherwise.

        :param corr_id: Correlation (RPC request) ID
        :type corr_id: int
        :returns: JSON-RPC2 compliant response from RPC server
        :rtype: dict
        """
        assert corr_id
        response_string = self.queue[corr_id]
        if response_string is None:
            return None
        response = {}
        try:
            response = json.loads(response_string)
        except Exception as error:
            raise InvalidJsonResponse(error)
        try:
            jsonrpc_utils.is_valid_response(response)
        except Exception as error:
            response = {
                "error": 400,
                "message": error,
                "data": response_string,
            }

        # Discard unused "id" field and return to user.
        # (we use rabbitmq correlation_id instead)
        if "id" in response and response["id"] is None:
            del response["id"]
        return response

    def request(self, method, params={}):
        """Format and send a JSON-RPC request, wait for response, and return result (which may indicate an error).

        :param method: The RPC method to call; this is not validated by the rpc client.
        :type method: str
        :param params: Any associated parameters to accompany the request (varies by method).
        :type params: dict
        :returns: JSON-RPC2 compliant response from RPC server; it may indicate error.
        :rtype: dict
        """
        # Send the request and store the requests' ID
        corr_id = self.send_request(method, params)

        # Wait for up to a maximum amount of time for a response.
        waited = 0
        response = {}
        while self.queue[corr_id] is None:
            waited = waited + self.wait_interval
            if waited > self.max_wait:
                # timeout error
                response["error"] = {"code": 500, "message": "Server timeout"}
                return response
            sleep(self.wait_interval)

        # Return the JSON-RPC2 response to the user (may be error).
        return self.get_response(corr_id)


class InvalidJsonResponse(Exception):
    pass


class ConfigurationError(Exception):
    pass


class ConnectionError(Exception):
    pass
