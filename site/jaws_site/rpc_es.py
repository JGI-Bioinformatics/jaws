import json
import logging

from amqpstorm import Message
from jaws_rpc import rpc_client


class RPCRequest(rpc_client.RpcClient):
    """Asynchronous remote procedure call (RPC) client class. This class inherits from the rpc_client.RpcClient
    module and have overridding methods for sending a json payload to a RMQ queue."""

    def __init__(self, rpc_params: dict, logger: logging) -> None:
        """Calls the parent rpc_client constructor to create a RMQ connection. The RMQ connection entries are passed
        into the rpc_entries parameter as a dictionary with the following required keys:
        {
            user: RabbitMQ user
            password: RabbitMQ password
            host: RabbitMQ host
            port: RabbitMQ port
            vhost: vhost name
            queue: queue name
        }

        Input paramters:
        :param rpc_params: dictionary containing the RMQ connections
        :type rpc_params: dict
        :param logger: log object
        :type logger: logging
        """
        super().__init__(rpc_params, logger)

    def send_request(self, payload: dict) -> str:
        """Format and send a JSON-RPC2 request but do not wait for a response.
        This function overwrites the rpc_cliennt.send_request method to allow a json payload
        to be passed in and used instead of passing in the method and params parameters.

        :param payload: The json payload to send to the RMQ queue.
        :type payload: dict
        :return: the ID of the RPC request (correlation ID)
        :rtype: int
        """
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

    def request(self, payload: dict) -> str:
        """Format and send a JSON-RPC request, wait for response, and return result (which may indicate an error).
        This function overwrites the rpc_cliennt.send_request method to allow a json payload
        to be passed in and used instead of passing in the method and params parameters.

        :param payload: The json payload to send to the RMQ queue.
        :type payload: dict
        :returns: JSON-RPC2 compliant response from RPC server; it may indicate error.
        :rtype: dict
        """
        # Convert dictionary to string
        payload = json.dumps(payload, default=str)

        # Send the request and store the requests' ID
        corr_id = self.send_request(payload)

        # Return the JSON-RPC2 response to the user (may be error).
        return self.get_response(corr_id)
