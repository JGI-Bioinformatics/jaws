"""
File contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
import os
import shutil


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-run.ini"
    content = """[LOCAL_RPC_SERVER]
host = localhost
vhost = jaws_test
queue = site_rpc
user = jaws
password = passw0rd1
num_threads = 5
max_retries = 3

[CENTRAL_RPC_SERVER]
host = currenthost
vhost = jaws_test
queue = eagle
user = jaws_eagle
password = succotash
num_threads = 5
max_retries = 3

[CENTRAL_RPC_CLIENT]
host = currenthost
vhost = jaws_test
queue = central_rpc
user = jaws_eagle
password = succotash

[GLOBUS]
client_id = foghorn_leghorn
endpoint_id = rooster
root_dir = cwd
default_dir = /

[DB]
dialect = mysql+mysqlconnector
host = myhost
port = 60032
user = elmer_fudd
password = hunting
db = hunting_sites

[CROMWELL]
url = http://localhost:8000

[SITE]
id = eagle
uploads_subdirectory = uploads
downloads_subdirectory = downloads
"""

    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture()
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-run.ini"
    content = """[CENTRAL_RPC_SERVER]
host = https://rmq.nersc.gov
vhost = jaws_test
user = bugs_bunny
password = xqweasdasa
max_retries = 10

[CENTRAL_RPC_CLIENT]
host = https://rmq.nersc.gov
user = bugs_bunny
password = xqweasdasa
vhost = jaws_test

[GLOBUS]
client_id = foghorn_leghorn
endpoint_id = rooster
root_dir = cwd
default_dir = /

[DB]
dialect = mysql+mysqlconnector
host = myhost
port = 60032
user = elmer_fudd
password = hunting
db = hunting_sites

[LOCAL_RPC_SERVER]
vhost = jaws_test

[CROMWELL]
url = http://localhost:8000

[SITE]
id = eagle
uploads_subdirectory = uploads
downloads_subdirectory = downloads
    """
    cfg.write_text(content)
    return cfg.as_posix()


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

    def ok(self):
        return False if self.status_code >= 400 else True

    def raise_for_status(self):
        if not self.ok:
            raise
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
                    "split_files": ["/home/jdoe/temp_aa", "/home/jdoe/temp_ad"]
                },
                "inputs": {"input_file": "/home/jdoe/cromwell/11.txt"},
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0",
                },
                "returnCode": 0,
                "backend": "Local",
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:55.000-05:00",
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
                    "continueOnReturnCode": "0",
                },
                "inputs": {"input_file": "f"},
                "backend": "Local",
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/0/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00",
            },
            {
                "executionStatus": "Done",
                "stdout": "/home/jdoe/shard-0/attempt-2/stdout",
                "shardIndex": 0,
                "outputs": {"count_file": "/home/jdoe/0/attempt-2/output.txt"},
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0",
                },
                "inputs": {"input_file": "f"},
                "returnCode": 0,
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/shard-0/attempt-2/stderr",
                "attempt": 2,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00",
            },
            {
                "executionStatus": "Done",
                "stdout": "/home/jdoe/shard-1/stdout",
                "shardIndex": 1,
                "outputs": {"count_file": "/home/jdoe/shard-1/output.txt"},
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0",
                },
                "inputs": {"input_file": "f"},
                "returnCode": 0,
                "backend": "Local",
                "end": "2016-02-04T13:47:56.000-05:00",
                "stderr": "/home/jdoe/shard-1/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00",
            },
        ],
        "sc_test.do_gather": [
            {
                "executionStatus": "Done",
                "stdout": "/home/jdoe/call-do_gather/stdout",
                "shardIndex": -1,
                "outputs": {"sum": 12},
                "runtimeAttributes": {
                    "failOnStderr": "true",
                    "continueOnReturnCode": "0",
                },
                "inputs": {
                    "input_files": [
                        "/home/jdoe/shard-0/attempt-2/output.txt",
                        "/home/jdoe/shard-1/output.txt",
                    ]
                },
                "returnCode": 0,
                "backend": "Local",
                "end": "2016-02-04T13:47:57.000-05:00",
                "stderr": "/home/jdoe/call-do_gather/stderr",
                "attempt": 1,
                "executionEvents": [],
                "start": "2016-02-04T13:47:56.000-05:00",
            }
        ],
    },
    "outputs": {
        "sc_test.do_gather.sum": 12,
        "sc_test.do_prepare.split_files": [
            "/home/jdoe/call-do_prepare/temp_aa",
            "/home/jdoe/call-do_prepare/temp_ad",
        ],
        "sc_test.do_scatter.count_file": [
            "/home/jdoe/shard-0/attempt-2/output.txt",
            "/home/jdoe/shard-1/output.txt",
        ],
    },
    "id": "8e592ed8-ebe5-4be0-8dcb-4073a41fe180",
    "inputs": {"sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"},
    "submission": "2016-02-04T13:47:55.000-05:00",
    "status": "Succeeded",
    "end": "2016-02-04T13:47:57.000-05:00",
    "start": "2016-02-04T13:47:55.000-05:00",
}

LOGS = {
    "id": "b3e45584-9450-4e73-9523-fc3ccf749848",
    "logs": {
        "call.ps": [
            {
                "stderr": "/home/user/call-ps/stderr6126967977036995110.tmp",
                "stdout": "/home/user/call-ps/stdout6128485235785447571.tmp",
            }
        ],
        "call.cgrep": [
            {
                "stderr": "/home/user/call-cgrep/stderr6126967977036995110.tmp",
                "stdout": "/home/user/call-cgrep/stdout6128485235785447571.tmp",
            }
        ],
        "call.wc": [
            {
                "stderr": "/home/user/call-wc/stderr6126967977036995110.tmp",
                "stdout": "/home/user/call-wc/stdout6128485235785447571.tmp",
            }
        ],
    },
}

ABORT = {"id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1", "status": "Aborted"}


@pytest.fixture()
def server_status_down():
    def get(url):
        return MockResponses(
            {"serviceName": "cori", "messages": "server unhealthy", "ok": "false"}, 500
        )

    return get


@pytest.fixture()
def server_status_up():
    def get(url):
        return MockResponses(
            {"serviceName": "cori", "messages": "server healthy", "ok": "true"}, 200
        )

    return get


@pytest.fixture()
def server_status_get():
    def get(url):
        return MockResponses({"id": CROMWELL_ID, "status": "Submitted"}, 200)

    return get


@pytest.fixture()
def log_file(tmp_path):
    logfile = tmp_path / "stderr.submit"
    content = ""
    for i in range(3000):
        content += f"this is line number {i}\n"
    logfile.write_text(content)

    return logfile.as_posix()


@pytest.fixture()
def user_dir(tmp_path):
    fpath = tmp_path / "user_dir"
    return fpath.as_posix()


@pytest.fixture()
def cromwell_run_dir(tmp_path):
    cromwell_path = tmp_path / "cromwell-execution"
    subdirs = [
        "call-asm_1",
        "call-asm_2",
        "call-asm_3",
        "call-callGenes",
        "call-circularizeAssembly",
        "call-createNonMitoReads",
        "call-doHmmSearch",
        "call-filterHighGc",
        "call-getContigsMatchingHmmSearch",
    ]

    # These dirs will have an rc file with non-zero exit code
    non_zero_dirs = ["call-asm_1", "call-circularizeAssembly", "call-filterHighGc"]

    # These dirs will have nonzero stderr.submit files (submission errors)
    submission_error_dirs = ["call-asm_1", "call-filterHighGc", "call-doHmmSearch"]

    for d in subdirs:
        dir_path = cromwell_path / d / "execution"
        dir_path.mkdir(parents=True)
        rc_file = dir_path / "rc"
        stdout = dir_path / "stdout"
        stderr = dir_path / "stderr"
        stdout_sub = dir_path / "stdout.submit"
        stderr_sub = dir_path / "stderr.submit"

        stdout.write_text(f"This is standard output from {d}")
        stdout_sub.write_text(f"This is submit stdout from {d}")
        if d in non_zero_dirs:
            rc_file.write_text("127")
            stderr.write_text(f"This is standard error. This {d} had an error")
            if d in submission_error_dirs:
                stderr_sub.write_text(f"This is submit stderr from {d}")
        else:
            rc_file.write_text("0")
            stderr.write_text(f"This is standard error. This {d} had no errors")
            stderr_sub.write_text("")
    return cromwell_path.as_posix()


@pytest.fixture()
def uploads_files():
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

    def _query_user_id(self, *args, **kwargs):
        return


class MockDb:
    def __init__(self):
        return

    def session(self):
        return MockSession()


class MockRun:
    """Mock Run object with useable defaults."""
    def __init__(self, **kwargs):
        self.user_id = kwargs.get("user_id", "jaws")
        self.upload_task_id = kwargs.get("upload_task_id", "1")
        self.submission_id = kwargs.get("submit_id", "2")
        self.cromwell_run_id = kwargs.get("cromwell_run_id", "myid")
        self.status = kwargs.get("status", "running")
        self.id = kwargs.get("id", "99")
        self.output_endpoint = kwargs.get("output_endpoint", "EXAMPLE_OUTPUT_ENDPOINT_ID")
        self.output_dir = kwargs.get("output_dir", ".")
        self.download_task_id = kwargs.get("download_task_id", "325")
        self.transfer_refresh_token = "EXAMPLE_GLOBUS_TRANSFER_TOKEN"
        self.email = "jaws@vog.gov"
        self.cromwell_workflow_dir = "/global/scratch/jaws/dev/cromwell-executions/test_wdl/myid"


class MockTransferClient:
    def __init__(self, status, transfer_result={"task_id": "325"}):
        self.status = status
        self.transfer_result = transfer_result

    def get_task(self, task_id):
        return self.status

    def submit_transfer(self, transfer_dat):
        return self.transfer_result


class MockTransferClientWithCopy:
    def __init__(self, status, transfer_result={"task_id": "325"}):
        self.status = status
        self.transfer_result = transfer_result

    def get_task(self, task_id):
        return self.status

    def submit_transfer(self, transfer_dat):
        print(f"copying {transfer_dat.src} to {transfer_dat.dst}")
        shutil.copytree(transfer_dat.src, transfer_dat.dst)
        return self.transfer_result


class MockUserQuery:
    def __init__(self, token):
        self.transfer_refresh_token = token


class MockUser:
    def __init__(self):
        self.id = "jaws_user"
        self.transfer_refresh_token = "1234567890"


def query_jaws_id(jawsd, run):
    return MockUser()


class MockTransferData:
    def __init__(self, *args, **kwargs):
        pass

    def add_item(self, src, dst, **kwargs):
        return


class MockTransferDataWithCopy:
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def add_item(self, src, dst, **kwargs):
        self.src = src
        self.dst = dst


@pytest.fixture()
def mock_query_user_id(monkeypatch):
    from jaws_run.daemon import Daemon

    monkeypatch.setattr(Daemon, "_query_user_id", query_jaws_id)
