"""
File contains all the mock classes and fixtures that will be used during testing.
"""
import pytest
import amqpstorm
import threading
import json


@pytest.fixture()
def amqp_config():
    return {
        "host": "localhost",
        "vhost": "vlocalhost",
        "user": "mrcromwell",
        "password": "secrets",
        "queue": "myqueue",
    }


@pytest.fixture()
def rpc_config():
    return {"num_threads": 5, "max_retries": 5}


@pytest.fixture()
def mock_connection(monkeypatch):
    """Do not attempt to create a connection to an AMQP server"""
    monkeypatch.setattr(amqpstorm, "UriConnection", MockConnection)


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
        # We want to be able to go through at least one iteration of starting our server
        # so we first set our flag as False to enter the loop. Then we turn it back to True
        # on the second iteration. It's a very ugly hack.
        if self.setting is True:
            self.setting = False
        elif self.setting is False:
            self.setting = True
        return self.setting


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
    """Mocks the RabbitMQ connection """

    def __init__(self, uri):
        self.uri = uri

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


def mock_status_get(url):
    cromwell_id = "15774623-0f76-49ef-828c-3aa0ccd024f5"
    return MockResponses({"id": cromwell_id, "status": "Submitted"}, 201)


def mock_metadata_get(url):
    cromwell_id = "15774623-0f76-49ef-828c-3aa0ccd024f5"
    return MockResponses(
        {
            "id": cromwell_id,
            "calls": {
                "backend": "Slurm",
                "backendLogs": {"logs": "path/to/logs"},
                "end": "123456789",
                "executionStatus": "Running",
                "start": "123456789",
            },
            "inputs": {"input": "align.wdl"},
            "outputs": {"output": "out"},
            "start": "1234556677890",
            "status": "Pending",
        },
        201,
    )


def mock_abort_get(url):
    cromwell_id = "15774623-0f76-49ef-828c-3aa0ccd024f5"
    return MockResponses({"id": cromwell_id, "status": "Aborted"}, 201)


def mock_labels_get(url):
    cromwell_id = "15774623-0f76-49ef-828c-3aa0ccd024f5"
    return MockResponses(
        {"labels": {"id": cromwell_id, "cromwell-workflow-id": "label-id-1"}}, 200
    )
