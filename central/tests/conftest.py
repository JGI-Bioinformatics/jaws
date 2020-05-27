"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
import amqpstorm
import threading


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-central.ini"
    content = """[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3306
user = jaws
password = passw0rd1
db = jaws

[GLOBUS]
client_id = ZZZZ

[SITE:LBNL]
amqp_host = rmq.foo.com
amqp_user = jaws
amqp_password = passw0rd2
amqp_vhost =
amqp_queue = jaws_rpc
globus_endpoint = XXXX
globus_basepath = "/global/scratch/jaws"
staging_subdir = "staging"
max_ram_gb = 1024

[SITE:NERSC]
amqp_host = rmq.bar.com
amqp_user = jaws
amqp_password = passw0rd3
amqp_vhost =
amqp_queue = jaws_rpc
globus_endpoint = YYYY
globus_basepath = "/"
staging_subdir = "/global/scratch/jaws/staging"
max_ram_gb = 2048
"""
    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture()
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
    content = """[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3305
user = j4ws
password = p455w0rd1
db = jaws

[GLOBUS]
client_id = ZZZZ
"""
    cfg.write_text(content)
    return cfg.as_posix()


rpc_client_dict = {
    "amqp_host": "rmq.foo.com",
    "amqp_user": "jaws",
    "amqp_password": "passw0rd2",
    "amqp_vhost": "",
    "amqp_queue": "jaws_rpc",
    "globus_endpoint": "XXXX",
    "globus_basepath": "\"/global/scratch/jaws\"",
    "staging_subdir": "staging",
    "max_ram_gb": 1024
}


@pytest.fixture()
def rpc_dict():
    return rpc_client_dict


@pytest.fixture()
def mock_connection(monkeypatch):
    """Do not attempt to create a connection to an AMQP server"""
    monkeypatch.setattr(amqpstorm, "Connection", MockConnection)


@pytest.fixture()
def mock_uri_connection(monkeypatch):
    """Do not attempt to create a connection to an AMQP server"""
    monkeypatch.setattr(amqpstorm, "UriConnection", MockUriConnection)


@pytest.fixture()
def mock_thread(monkeypatch):
    monkeypatch.setattr(threading, "Thread", MockThread)


@pytest.fixture()
def mock_message(monkeypatch):
    monkeypatch.setattr(amqpstorm, "Message", MockMessage)


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

    def declare(self, queue_name='', exclusive=False):
        return {"jsonrpc": "2.0", "result": "result", "error": "error",
                "id": "id", "queue": queue_name}


class MockBasic:
    def __init__(self):
        pass

    def consume(self, consumer, no_ack=False, queue=''):
        return "consumer tag"

    def qos(self):
        return {"jsonrpc": "2.0", "result": "result", "id": "id"}

    def publish(self, body, routing_key, exchange, properties, mandatory, immediate):
        return True


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

    def __init__(self, amqp_host, amqp_user, amqp_password):
        self.host = amqp_host
        self.user = amqp_user
        self.password = amqp_password

    def is_closed(self):
        return False

    def check_for_errors(self):
        return

    def channel(self):
        return MockChannel()


class MockUriConnection:
    """Mocks the RabbitMQ connection """

    def __init__(self, uri):
        self.uri = uri

    def is_closed(self):
        return False

    def check_for_errors(self):
        return

    def channel(self):
        return MockChannel()


class MockErrorConnection:
    def __init__(self, uri):
        raise amqpstorm.AMQPConnectionError(f"Could not connect to {uri}")


class MockMessage:
    def __init__(self, channel, payload):
        self.correlation_id = "corr_id"
        self.body = payload

    def ack(self):
        return

    def channel(self):
        return MockChannel()

    @staticmethod
    def create(channel, body, properties=None):
        return MockMessage(channel, body)

    def publish(self, routing_key):
        return
