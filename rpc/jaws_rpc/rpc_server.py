"""
A Scalable and threaded Consumer that will automatically re-connect on failure.
"""
import threading
import time
import json
import amqpstorm
from sqlalchemy.orm import scoped_session


DEFAULT_PORT = 5672
DEFAULT_NUM_THREADS = 5
DEFAULT_MAX_RETRIES = 3


class InvalidRequest(Exception):
    pass


class InvalidResponse(Exception):
    pass


class Consumer(object):
    def __init__(self, logger, queue, operations, Session=None):
        """Initialize Consumer object

        :param queue: The name of the queue from which to retrieve messages.
        :type queue: str
        :return:
        """
        self.queue = queue
        self.operations = operations
        self.Session = Session
        self.logger = logger
        self.channel = None
        self.active = False

    def start(self, connection):
        """Start the consumer"""
        self.channel = None
        try:
            self.active = True
            self.channel = connection.channel(rpc_timeout=10)
            self.channel.basic.qos(1)
            self.channel.queue.declare(self.queue)
            self.channel.basic.consume(self, self.queue, no_ack=False)
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
        # DECODE AND VALIDATE REQUEST
        try:
            request = self.validate_request(message.body)
        except InvalidRequest as error:
            self.logger.error(f"Invalid JSON-RPC2 request; {error}: {request}")
            response = {
                "jsonrpc": "2.0",
                "error": {"code": 400, "message": f"Invalid RPC request: {request}"},
            }
            return self.__respond__(message, response)

        # GET METHOD AND OPTIONAL PARAMS
        method = request["method"]
        params = request["params"] if "params" in request else {}

        # VALIDATE METHOD AND PARAMS
        if method not in self.operations:
            self.logger.error(f"Unknown JSON-RPC2 method: {method}")
            response = {
                "jsonrpc": "2.0",
                "error": {"code": 400, "message": f"Unknown RPC method: {method}"},
            }
            return self.__respond__(message, response)
        if "required_params" in self.operations:
            for required_param in self.operations["required_params"]:
                if required_param not in params:
                    self.logger.error(f"Method {method} missing {required_param}")
                    response = {
                        "jsonrpc": "2.0",
                        "error": {
                            "code": 400,
                            "message": f"Method {method} missing {required_param}",
                        },
                    }
                    return self.__respond__(message, response)

        # GET RESPONSE FROM DISPATCH FUNCTION
        self.logger.debug(f"RPC method {method} with {params}")
        proc = self.operations[method]["function"]
        # rpc procedures are either called with a db session or not, depending on whether the
        # rpc manager was initialized with a scoped_session obj or not, because an rpc server
        # may or may not have an associated (sqlalchemy) db.
        response = None
        session = None
        if self.Session:
            session = self.Session()
            response = proc(params, session)
        else:
            response = proc(params)
        self.__respond__(message, response)
        if session:
            Session.remove()

    def __respond__(self, message, response):
        try:
            self.validate_response(response)
        except InvalidResponse as error:
            self.logger.error(f"Invalid response; {error}: {response}")
            response = {
                "jsonrpc": "2.0",
                "error": {
                    "code": 500,
                    "message": f"Method returned invalid response; {error}: {response}",
                },
            }

        # INIT RESPONSE MESSAGE OBJ
        properties = {"correlation_id": message.correlation_id}
        response_message = amqpstorm.Message.create(
            message.channel, json.dumps(response), properties
        )

        # DELIVER RESPONSE AND REMOVE REQUEST MESSAGE FROM QUEUE
        response_message.publish(message.reply_to)
        message.ack()

    def validate_request(self, message_str):
        """Verifies the request dict conforms to JSON-RPC2 spec.  Raises error if invalid.

        :param request: The message body with json-encoded request
        :type request: str
        :return: valid json-rpc2 request dict; raise otherwise
        :rtype: dict
        """
        try:
            request = json.loads(message_str)
        except Exception:
            raise InvalidRequest("not valid JSON")
        if request is None:
            raise InvalidRequest("request is undefined")
        elif type(request) is not dict:
            raise InvalidRequest("request is not a dict")
        elif "method" not in request:
            raise InvalidRequest("method is undefined")
        elif "jsonrpc" not in request:
            raise InvalidRequest("jsonrpc is undefined")
        elif request["jsonrpc"] != "2.0":
            raise InvalidRequest("jsonrpc is not 2.0")
        else:
            request_valid_keys = ["jsonrpc", "id", "method", "params"]
            for key in request.keys():
                if key in request_valid_keys:
                    pass
                else:
                    raise InvalidRequest(f"invalid key, {key}")
        return request

    def validate_response(self, response):
        """Verifies the response dict conforms to JSON-RPC2 spec.  Raises error if invalid.

        :param response: The response object to validate
        :type response: dict
        :return: True if valid; raise otherwise
        :rtype: bool
        """
        if "jsonrpc" not in response:
            raise InvalidResponse("jsonrpc is undefined")
        elif response["jsonrpc"] != "2.0":
            raise InvalidResponse("jsonrpc is not 2.0")
        elif "error" in response:
            if type(response["error"]) is not dict:
                raise InvalidResponse("error is not a dict")
            elif "code" not in response["error"]:
                raise InvalidResponse("error code is undefined")
            elif "message" not in response["error"]:
                raise InvalidResponse("error message is undefined")
            elif "result" in response:
                raise InvalidResponse("cannot have error and result")
            else:
                response_error_valid_keys = ["code", "message", "data"]
                for key in response["error"].keys():
                    if key in response_error_valid_keys:
                        pass
                    else:
                        raise InvalidResponse(f"invalid error key, {key}")
        elif "result" not in response:
            raise InvalidResponse("neither result or code is defined")
        else:
            response_valid_keys = ["jsonrpc", "id", "error", "result"]
            for key in response.keys():
                if key in response_valid_keys:
                    pass
                else:
                    raise InvalidResponse(f"invalid key, {key}")
        return True


class RpcServer(object):
    def __init__(self, params, logger, operations, Session=None) -> None:
        """
        Init Rpc Server
        :param params: configuration parameters
        :ptype params: dict
        :param logger: logging object
        :ptype logger: Logging
        :param operations: valid operations required params and dispatch table
        :ptype operations: dict
        :param Session: thread-local session object factory
        :ptype Session: sqlalchemy.orm.scoped_session
        """
        self.logger = logger
        self.params = {}
        for required_param in ["host", "vhost", "user", "password", "queue"]:
            if required_param not in params:
                raise ConfigurationError(f"{required_param} required")
            self.params[required_param] = params[required_param]
        self.num_threads = int(params.get("num_threads", DEFAULT_NUM_THREADS))
        self.max_retries = int(params.get("max_retries", DEFAULT_MAX_RETRIES))
        self.params["port"] = int(params.get("port", DEFAULT_PORT))
        self.logger.info(
            f"Connecting to host:{params['host']}, vhost:{params['vhost']}, queue:{params['queue']}"
        )
        self.operations = operations
        self.Session = Session
        self.consumers = [
            Consumer(self.logger, params["queue"], self.operations, self.Session)
            for _ in range(self.num_threads)
        ]
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
            consumer = Consumer(
                self.params["queue"], self.operations, self.Session
            )
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
                return amqpstorm.Connection(
                    self.params["host"],
                    self.params["user"],
                    self.params["password"],
                    int(self.params["port"]),
                    virtual_host=self.params["vhost"],
                    heartbeat=10,
                )
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


class ConfigurationError(Exception):
    pass
