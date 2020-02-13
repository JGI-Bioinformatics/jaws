"""
AmpqStorm RPC Client for interacting with JAWS-Sites for analysis-run management.

Sites may use either the same or different AMPQ servers, as specified in their config files.
"""

import threading
import amqpstorm
from amqpstorm import Message
import json
from time import sleep
import urllib.parse
#from flask import current_app

import config

class RPC_Client(object):
    """Asynchronous Rpc client."""

    def __init__(self, host, user, password, queue, vhost=None):
        self.queue = {}
        self.host = host
        self.vhost = vhost
        self.username = user
        self.password = password
        self.channel = None
        self.connection = None
        self.callback_queue = None
        self.rpc_queue = queue
        self.open()

    def open(self):
        """Open Connection."""
        if self.vhost:
            uri = "amqp://" + self.username + ":" + urllib.parse.quote_plus(self.password) + "@" + self.host + ":5672/" + self.vhost + "?heartbeat=60"
            self.connection = amqpstorm.UriConnection(uri)
        else:
            self.connection = amqpstorm.Connection(self.host, self.username, self.password)
        self.channel = self.connection.channel()
        self.channel.queue.declare(self.rpc_queue)
        result = self.channel.queue.declare(exclusive=True)
        self.callback_queue = result['queue']
        self.channel.basic.consume(self._on_response, no_ack=True, queue=self.callback_queue)
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
        """
        Format and send a JSON-RPC request and return the request's correlation_id, but do not wait for a response.
        """
        # Construct JSON-RPC2 payload string
        query = {
            "jsonrpc" : "2.0",
            "method" : method,
            "params" : params
        }
        payload = json.dumps(query)

        # Create the Message object.
        message = Message.create(self.channel, payload)
        message.reply_to = self.callback_queue

        # Create an entry in our local dictionary, using the automatically
        # generated correlation_id as our key.
        self.queue[message.correlation_id] = None

        # Publish the RPC request.
        message.publish(routing_key=self.rpc_queue)

        # Return the Unique ID used to identify the request.
        return message.correlation_id


    def get_response(self, corr_id):
        """
        Return the JSON-RPC2 response if exists, None otherwise.
        """
        assert(corr_id)
        response_string = self.queue[corr_id]
        if response_string is None: return None
        response = {}
        try:
            response = json.loads(response_string)
        except:
            response = {}
            response["error"] = {
                "code" : 500,
                "message" : "Invalid response: %s" % (response_string,)
            }
        return response


    def request(self, method, params={}, max_wait = 10):
        """
        Format and send a JSON-RPC request, wait for response, and return result (which may indicate an error).
        """
        # Send the request and store the requests' ID
        corr_id = self.send_request(method, params)

        # Wait for up to a maximum amount of time for a response.
        wait_interval = 0.25 # seconds
        waited = 0
        response = {}
        while self.queue[corr_id] is None:
            waited = waited + wait_interval    
            if waited > max_wait:
                # timeout error
                response["error"] = {
                    "code" : 500,
                    "message": "Server timeout"
                }
                return response
            sleep(wait_interval)

        # Return the response to the user (may be error).
        return self.get_response(corr_id)

