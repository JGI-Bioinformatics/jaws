"""
File contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
import amqpstorm
import threading
import json
import os
import shutil


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
    content = """[AMQP]
host = currenthost
vhost = current_vhost
user = daffy_duck
password = succotash
queue = rabbit_season

[RPC]
num_threads = 5
max_retries = 3

[GLOBUS]
client_id = foghorn_leghorn
endpoint_id = rooster
root_dir = cwd

[DB]
dialect = mysql+mysqlconnector
host = myhost
port = 60032
user = elmer_fudd
password = hunting
db = hunting_sites

[CROMWELL]
workflows_url = http://localhost:8000/api/workflows/v1
engine_status_url = http://localhost:8000/engine/v1/status

[SITE]
id = eagle
staging_subdirectory = staging
results_subdirectory = results
"""

    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture()
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
    content = """[AMQP]
host = https://rmq.nersc.gov
user = bugs_bunney
password = xqweasdasa

[RPC]
max_retries = 10

[CROMWELL]
workflows_url = http://localhost:8000/api/workflows/v1
engine_status_url = http://localhost:8000/engine/v1/status
    """
    cfg.write_text(content)
    return cfg.as_posix()


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
        # We want to be able to go through at least one iteration of
        # starting our server so we first set our flag as False to enter the
        # loop. Then we turn it back to True
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
        return {"jsonrpc": "2.0", "result": "result", "error": "error",
                "id": "id"}


class MockBasic:
    def __init__(self):
        pass

    def consume(self, consumer, rpc, no_ack=False):
        return "consumer tag"

    def qos(self, number):
        return {"jsonrpc": "2.0", "result": "result", "error": "error",
                "id": "id"}


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


CROMWELL_ID = "15774623-0f76-49ef-828c-3aa0ccd024f5"

METADATA = {
    "workflowName": "sc_test",
    "calls": {
        "sc_test.do_prepare": [
            {
                "executionStatus": "Done",
                "stdout": "/path/to/stdout",
                "shardIndex": -1,
                "outputs": {
                    "split_files": [
                        "/home/jdoe/temp_aa",
                        "/home/jdoe/temp_ad"
                    ]
                },
                "inputs": {
                    "input_file": "/home/jdoe/cromwell/11.txt"
                },
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0"
                },
                "returnCode": 0,
                "backend": "Local",
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:55.000-05:00"
            }
        ],
        "sc_test.do_scatter": [
            {
                "executionStatus": "Preempted",
                "stdout": "/home/shard-0/stdout",
                "shardIndex": 0,
                "outputs": {},
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0"
                },
                "inputs": {
                    "input_file": "f"
                },
                "backend": "Local",
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/0/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00"
            },
            {
                "executionStatus": "Done",
                "stdout": "/home/jdoe/shard-0/attempt-2/stdout",
                "shardIndex": 0,
                "outputs": {
                    "count_file": "/home/jdoe/0/attempt-2/output.txt"
                },
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0"
                },
                "inputs": {
                    "input_file": "f"
                },
                "returnCode": 0,
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/shard-0/attempt-2/stderr",
                "attempt": 2,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00"
            },
            {
                "executionStatus": "Done",
                "stdout": "/home/jdoe/shard-1/stdout",
                "shardIndex": 1,
                "outputs": {
                    "count_file": "/home/jdoe/shard-1/output.txt"
                },
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0"
                },
                "inputs": {
                    "input_file": "f"
                },
                "returnCode": 0,
                "backend": "Local",
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/shard-1/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00"
            }
        ],
        "sc_test.do_gather": [
            {
                "executionStatus": "Done",
                "stdout": "/home/jdoe/call-do_gather/stdout",
                "shardIndex": -1,
                "outputs": {
                    "sum": 12
                },
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0"
                },
                "inputs": {
                    "input_files": [
                        "/home/jdoe/shard-0/attempt-2/output.txt",
                        "/home/jdoe/shard-1/output.txt"
                    ]
                },
                "returnCode": 0,
                "backend": "Local",
                "end": "2016-02-04T13:47:57.000-05:00",
                "stderr": "/home/jdoe/call-do_gather/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00"
            }
        ]
    }, "outputs": {
        "sc_test.do_gather.sum": 12,
        "sc_test.do_prepare.split_files": [
            "/home/jdoe/call-do_prepare/temp_aa",
            "/home/jdoe/call-do_prepare/temp_ad"
        ],
        "sc_test.do_scatter.count_file": [
            "/home/jdoe/shard-0/attempt-2/output.txt",
            "/home/jdoe/shard-1/output.txt"
        ]}, "id": "8e592ed8-ebe5-4be0-8dcb-4073a41fe180", "inputs": {
        "sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"},
    "submission": "2016-02-04T13:47:55.000-05:00", "status": "Succeeded",
    "end": "2016-02-04T13:47:57.000-05:00", "start":
        "2016-02-04T13:47:55.000-05:00"
}

LOGS = {
    "id": "b3e45584-9450-4e73-9523-fc3ccf749848",
    "logs": {
        "call.ps": [
            {
                "stderr": "/home/user/call-ps/stderr6126967977036995110.tmp",
                "stdout": "/home/user/call-ps/stdout6128485235785447571.tmp"
            }
        ],
        "call.cgrep": [
            {
                "stderr": "/home/user/call-cgrep/stderr6126967977036995110.tmp",
                "stdout": "/home/user/call-cgrep/stdout6128485235785447571.tmp"
            }
        ],
        "call.wc": [
            {
                "stderr": "/home/user/call-wc/stderr6126967977036995110.tmp",
                "stdout": "/home/user/call-wc/stdout6128485235785447571.tmp"
            }
        ]
    }
}

ABORT = {
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "status": "Aborted"
}


@pytest.fixture()
def server_status_down():
    def get(url):
        return MockResponses({"serviceName": "cori",
                              "messages": "server unhealthy",
                              "ok": "false"}, 500)

    return get


@pytest.fixture()
def server_status_up():
    def get(url):
        return MockResponses({"serviceName": "cori",
                              "messages": "server healthy",
                              "ok": "true"}, 200)

    return get


@pytest.fixture()
def server_status_get():
    def get(url):
        return MockResponses({"id": CROMWELL_ID, "status": "Submitted"}, 200)

    return get


@pytest.fixture
def metadata_get():
    def get(url):
        return MockResponses(METADATA, 200)

    return get


@pytest.fixture()
def logs_get():
    def get(url):
        return MockResponses(LOGS, 200)

    return get


@pytest.fixture()
def abort_post():
    def post(url):
        return MockResponses({"id": CROMWELL_ID, "status": "Aborted"}, 200)

    return post


@pytest.fixture()
def log_file(tmp_path):
    logfile = tmp_path / "stderr.submit"
    content = ""
    for i in range(3000):
        content += f"this is line number {i}\n"
    logfile.write_text(content)

    return logfile.as_posix()


@pytest.fixture()
def cromwell_run_dir(tmp_path):
    subdirs = ["call-asm_1", "call-asm_2", "call-asm_3",
               "call-callGenes", "call-circularizeAssembly",
               "call-createNonMitoReads", "call-doHmmSearch",
               "call-filterHighGc", "call-getContigsMatchingHmmSearch"]

    # These dirs will have an rc file with non-zero exit code
    non_zero_dirs = ["call-asm_1", "call-circularizeAssembly",
                     "call-filterHighGc"]

    for d in subdirs:
        dir_path = tmp_path / d / "execution"
        dir_path.mkdir(parents=True)
        rc_file = dir_path / "rc"
        if d in non_zero_dirs:
            rc_file.write_text("127")
        else:
            rc_file.write_text("0")
        stdout = dir_path / "stdout"
        stderr = dir_path / "stderr"
        stdout_sub = dir_path / "stdout.submit"
        stderr_sub = dir_path / "stderr.submit"
        stdout.write_text(f"This is standard output from {d}")
        stderr.write_text(f"This is standard error. This {d} had an error")
        stdout_sub.write_text(f"This is submit stdout from {d}")
        stderr_sub.write_text(f"This is submit stderr from {d}")
    return tmp_path.as_posix()


@pytest.fixture()
def staging_files():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "1")

    os.mkdir(root_dir)

    for f in ["2.wdl", "2.json", "2.zip"]:
        file_path = os.path.join(root_dir, f)
        with open(file_path, "w") as outfile:
            outfile.write(f"output for {f}")

    yield

    shutil.rmtree(root_dir)


class MockSession:

    def __init__(self):
        return

    def commit(self):
        return

    def close_all(self):
        return


class MockDb:

    def __init__(self):
        return

    def session(self):
        return MockSession()


class MockRun:

    def __init__(self, user_id, task_id, submit_id, cromwell_id, id, status):
        self.user_id = user_id
        self.upload_task_id = task_id
        self.submission_id = submit_id
        self.cromwell_id = cromwell_id
        self.id = id
        self.status = status
        self.output_endpoint = "."
        self.download_task_id = "123"


class MockTransferClient:

    def __init__(self, status, transfer_result={"task_id": "325"}):
        self.status = status
        self.transfer_result = transfer_result

    def get_task(self, task_id):
        return self.status

    def submit_transfer(self, transfer_dat):
        return self.transfer_result


class MockUser:

    def __init__(self):
        self.id = "jaws_user"
        self.transfer_refresh_token = "1234567890"


def query_jaws_id(jawsd, run):
    return MockUser()


class MockTransferData:

    def __init__(self, *args, **kwargs):
        return

    def add_item(self, src, dest, **kwargs):
        return


@pytest.fixture()
def mock_query_user_id(monkeypatch):
    from jaws_site.jawsd import Daemon
    monkeypatch.setattr(Daemon, "_query_user_id", query_jaws_id)
