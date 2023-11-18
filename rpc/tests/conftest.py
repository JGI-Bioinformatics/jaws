import configparser
import json
import os
import threading

import amqpstorm
import pytest


@pytest.fixture
def configuration():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    config = configparser.Configparser()
    config.read(f"{base_dir}/rpc.ini")
    return config


@pytest.fixture()
def mock_connection(monkeypatch):
    """Do not attempt to create a connection to an AMQP server"""
    monkeypatch.setattr(amqpstorm, "Connection", MockConnection)


@pytest.fixture()
def mock_thread(monkeypatch):
    monkeypatch.setattr(threading, "Thread", MockThread)


@pytest.fixture()
def mock_event(monkeypatch):
    monkeypatch.setattr(threading, "Event", MockEvent)


@pytest.fixture()
def mock_message(monkeypatch):
    monkeypatch.setattr(amqpstorm, "Message", MockMessage)


class MockEvent:
    def __init__(self):
        self.setting = True

    def clear(self):
        return

    def is_set(self):
        # We want to be able to go through at least one iteration of
        # starting our server so we first set our flag as False to enter the
        # loop. Then we turn it back to True
        # on the second iteration. It's a very ugly hack.
        if self.setting is True:
            self.setting = False
        elif self.setting is False:
            self.setting = True
        return self.setting

    def set(self):
        return


class MockThread:
    """
    Create a mock threading class so we can control when to start and stop without
    having to wait for a threading event
    """

    def __init__(self, target=None, args=None):
        self.target = target
        self.args = args

    def start(self):
        return


class MockQueue:
    def __init__(self):
        pass

    def declare(self, queue_name):
        return {"jsonrpc": "2.0", "result": "result", "error": "error", "id": "id"}


class MockBasic:
    def __init__(self):
        pass

    def consume(self, consumer, rpc, no_ack=False):
        return "consumer tag"

    def qos(self, number):
        return {"jsonrpc": "2.0", "result": "result", "error": "error", "id": "id"}


class MockChannel:
    def __init__(self):
        self.basic = MockBasic()
        self.queue = MockQueue()

    def basic(self):
        return self.basic

    def queue(self):
        return self.queue

    def start_consuming(self):
        return

    def consumer_tags(self):
        return ["consumer_tag"]

    def close(self):
        return


class MockConnection:
    """Mocks the RabbitMQ connection"""

    def __init__(self, host, user, password, port, virtual_host="", heartbeat=0):
        self.host = host
        self.user = user
        self.password = password
        self.port = port
        self.vhost = virtual_host
        self.heartbeat = heartbeat

    def is_closed(self):
        return False

    def check_for_errors(self):
        return

    def channel(self, rpc_timeout):
        return MockChannel()


class MockErrorConnection:
    def __init__(self, uri):
        raise amqpstorm.AMQPConnectionError(f"Could not connect to {uri}")


class MockMessage:
    def __init__(self, method, params, corr_id):
        self.method = method
        self.params = params
        self.correlation_id = corr_id
        self.body = json.dumps({"method": method, "params": params})
        self.reply_to = "other Message"

    def ack(self):
        return

    def channel(self):
        return MockChannel()

    @staticmethod
    def create(channel, response, properties):
        return MockResponses(response, 200)

    def publish(self, reply_to):
        return


class MockResponses:
    def __init__(self, json_data, status_code):
        self.status_code = status_code
        self.json_data = json_data

    def status_code(self):
        return self.status_code

    def json(self):
        return self.json_data

    def publish(self, reply_to):
        return
