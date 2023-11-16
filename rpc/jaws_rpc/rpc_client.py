"""RpcClient sends JSON-RPC2 messages to an RpcServer."""

import threading
import amqpstorm
from amqpstorm import Message
import json
from time import sleep
from jaws_rpc import jsonrpc_utils


DEFAULT_PORT = 5672
DEFAULT_WAIT_INTERVAL = 0.25
DEFAULT_MAX_WAIT = 10
DEFAULT_MESSAGE_TTL = 10  # expires in seconds or 0=doesn't expire


class RpcClient(object):
    """Asynchronous remote procedure call (RPC) client class."""

    def __init__(self, params, logger):
        """Constructor

        :param params: A dictionary containing configuration parameters.
        :type params: dict
        """
        self.params = {}
        self.logger = logger
        self.queue = {}
        self.channel = None
        self.connection = None
        self.callback_queue = None
        for required_param in ["host", "vhost", "user", "password", "queue"]:
            if required_param not in params:
                raise ConfigurationError(f"{required_param} required")
            self.params[required_param] = params[required_param]
        self.params["port"] = int(params.get("port", DEFAULT_PORT))
        self.wait_interval = float(
            params.get("rpc_wait_interval", DEFAULT_WAIT_INTERVAL)
        )
        self.max_wait = int(params.get("rpc_max_wait", DEFAULT_MAX_WAIT))
        self.message_ttl = int(params.get("rpc_message_ttl", DEFAULT_MESSAGE_TTL))
        if self.message_ttl > self.max_wait:
            raise ConfigurationError("rpc_message_ttl must be <= rpc_max_wait")
        self.open()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close connection to RabbitMQ when exiting a context manager (with statement)"""
        if self.channel is not None:
            self.channel.close()
            self.channel = None
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def __del__(self):
        """Close connection to RabbitMQ when object is destroyed"""
        if self.channel is not None:
            self.channel.close()
        if self.connection is not None:
            self.connection.close()

    def open(self):
        """Open connection to RabbitMQ"""
        self.logger.debug(
            f"Open connection to {self.params['host']}:{self.params['queue']}"
        )
        try:
            self.connection = amqpstorm.Connection(
                self.params["host"],
                self.params["user"],
                self.params["password"],
                int(self.params["port"]),
                virtual_host=self.params["vhost"],
            )
        except Exception as error:
            raise ConnectionError(error)
        self.channel = self.connection.channel()
        self.channel.queue.declare(self.params["queue"])
        result = self.channel.queue.declare(
            durable=False, exclusive=False, auto_delete=True
        )
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
        properties = {}
        if self.message_ttl:
            properties["expiration"] = str(self.message_ttl)
        message = Message.create(self.channel, payload, properties=properties)
        message.reply_to = self.callback_queue

        # Create an entry in our local dictionary, using the automatically
        # generated correlation_id as our key.
        self.queue[message.correlation_id] = None

        # Publish the RPC request.
        self.logger.debug(
            f"Publishing message {message.correlation_id} to {self.params['queue']}"
        )
        try:
            message.publish(routing_key=self.params["queue"])
        except Exception as error:
            raise error

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
                response["error"] = {
                    "code": 408,
                    "message": "Server timeout: The service is unable to respond at this time; please try again later.",
                }
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
