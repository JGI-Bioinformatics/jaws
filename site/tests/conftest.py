"""
File contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
import os
import shutil
from pathlib import Path
from dataclasses import dataclass


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
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
client_id = AAAA
client_secret = BBBB
endpoint_id = rooster
host_path = /global/scratch/jaws
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
uploads_dir = /global/scratch/jaws/jaws-dev/uploads
[AWS]
aws_access_key_id = AAAA
aws_secret_access_key = BBBB
s3_bucket = CCCC
"""

    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
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
client_id = AAAA
client_secret = BBBB
endpoint_id = rooster
host_path = /global/scratch/jaws
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
uploads_dir = /global/scratch/jaws/jaws-dev/uploads
[AWS]
aws_access_key_id = AAAA
aws_secret_access_key = BBBB
s3_bucket = CCCC
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
    root_dir = os.path.join(home_dir, "XXXX")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    for f in ["XXXX.wdl", "XXXX.json", "XXXX.orig.json", "XXXX.zip"]:
        file_path = os.path.join(root_dir, f)
        with open(file_path, "w") as outfile:
            outfile.write(f"output for {f}")

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def uploads_files_without_zip():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "WWWW")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    for f in ["WWWW.wdl", "WWWW.json", "WWWW.orig.json"]:
        file_path = os.path.join(root_dir, f)
        with open(file_path, "w") as outfile:
            outfile.write(f"output for {f}")

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def uploads_files_missing_json():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "YYYY")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    file_path = os.path.join(root_dir, "YYYY.wdl")
    with open(file_path, "w") as outfile:
        outfile.write(f"workflow test { ... }")

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def uploads_files_empty_wdl():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "ZZZZ")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    file_path = os.path.join(root_dir, "ZZZZ.wdl")
    Path(file_path).touch()

    file_path = os.path.join(root_dir, "ZZZZ.json")
    with open(file_path, "w") as outfile:
        outfile.write("{}")

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def transfer_dirs(tmp_path):
    uploads_dir = tmp_path / "cwd/uploads/jaws"

    uploads_dir.mkdir(parents=True)

    for f in ["XXXX.wdl", "XXXX.json", "XXXX.orig.json", "XXXX.zip"]:
        file_path = uploads_dir / f
        file_path.write_text(f"output for {f}")

    yield

    shutil.rmtree(uploads_dir)


class MockSession:
    def __init__(self):
        return

    def commit(self):
        return

    def close(self):
        return

    def close_all(self):
        return

    def _query_user_id(self, *args, **kwargs):
        return


class MockRunModel:
    """Mock Run sqlalchemy orm model object with useable defaults."""

    def __init__(self, **kwargs):
        self.user_id = kwargs.get("user_id", "jaws")
        self.upload_task_id = kwargs.get("upload_task_id", "1")
        self.submission_id = kwargs.get("submission_id", "XXXX")
        self.cromwell_run_id = kwargs.get("cromwell_run_id", "myid")
        self._status = kwargs.get("status", "running")
        self.id = kwargs.get("id", "99")
        self.output_endpoint = kwargs.get(
            "output_endpoint", "EXAMPLE_OUTPUT_ENDPOINT_ID"
        )
        self.output_dir = kwargs.get("output_dir", ".")
        self.download_task_id = kwargs.get("download_task_id", "325")
        self.email = "jaws@vog.gov"
        self.cromwell_workflow_dir = (
            "/global/scratch/jaws/dev/cromwell-executions/test_wdl/myid"
        )

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, new_status: str):
        self._status = new_status


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


class MockUser:
    def __init__(self):
        self.id = "jaws_user"


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
    from jaws_site.daemon import Daemon

    monkeypatch.setattr(Daemon, "_query_user_id", query_jaws_id)


@pytest.fixture()
def mock_data_transfer(monkeypatch):
    # import jaws_site.datatansfer_protocol
    from jaws_site.runs import Run
    from jaws_site.datatransfer_protocol import (
        Status,
        DataTransferError,
        DataTransferAPIError,
        DataTransferNetworkError,
        DataTransferFactory,
    )

    class MockDataTransfer:
        @staticmethod
        def submit_upload(metadata: str, manifest_file: list):
            if data_obj.raises.DataTransferError:
                raise DataTransferError()
            elif data_obj.raises.DataTransferAPIError:
                raise DataTransferAPIError()
            elif data_obj.raises.DataTransferNetworkError:
                raise DataTransferNetworkError()
            return "123"

        @staticmethod
        def submit_download(metadata: dict, src_dir: str, dst_dir: str):
            if data_obj.raises.DataTransferError:
                raise DataTransferError()
            elif data_obj.raises.DataTransferAPIError:
                raise DataTransferAPIError()
            elif data_obj.raises.DataTransferNetworkError:
                raise DataTransferNetworkError()
            return "456"

        @staticmethod
        def transfer_status(task_id: str):
            if data_obj.raises.DataTransferError:
                raise DataTransferError()
            elif data_obj.raises.DataTransferAPIError:
                raise DataTransferAPIError()
            elif data_obj.raises.DataTransferNetworkError:
                raise DataTransferNetworkError()

            if data_obj.status['succeeded']:
                return Status.succeeded
            elif data_obj.status['failed']:
                return Status.failed
            elif data_obj.status['transferring']:
                return Status.transferring
            elif data_obj.status['inactive']:
                return Status.inactive

        @staticmethod
        def cancel_transfer(task_id: str):
            if data_obj.raises.DataTransferError:
                raise DataTransferError()
            elif data_obj.raises.DataTransferAPIError:
                raise DataTransferAPIError()
            elif data_obj.raises.DataTransferNetworkError:
                raise DataTransferNetworkError()

    def mock_data_transfer(*args, **kwargs):
        return MockDataTransfer()

    def mock_get_data_transfer_type(*args, **kwargs):
        return 'globus_transfer'

    monkeypatch.setattr(DataTransferFactory, '__new__', mock_data_transfer)
    monkeypatch.setattr(Run, '_get_data_transfer_type', mock_get_data_transfer_type)

    @dataclass
    class DataTransferExceptions():
        DataTransferError = False
        DataTransferAPIError = False
        DataTransferNetworkError = False

    class Data():
        raises = DataTransferExceptions()
        status = {
            'succeeded': False,
            'failed': False,
            'transferring': False,
            'inactive': False,
        }

    data_obj = Data()
    return data_obj


@pytest.fixture
def mock_db_session():
    """Fixture to mockup the sqlalchemy session.

    session = database.session()  # here, session is a sqlalchemy sessionmaker object.
    result = session.query(table).\
        filter_by(**filters).\
        order_by(table.start_date.desc()).\
        all()

    will be mocked up. The result variable can be set in the pytest function using the mock_db.query variable.
    This variable accepts a list of either dictionaries, or another list of dictionaries.

    Ex: if the test function performs one sqlalchemy query, and we want to specify the return value of that query,
    (i.e., one query returning two entries):
    mock_db_session.output(
        [
            {'employee': 'John', 'Title': 'analyst'},
            {'employee': 'Mary', 'Title': 'PI'},
        ]
    )

    If the test function performs multiple sqlalchemy queries and we want to return the same result, add
    repeat=True to the mock_db_session.output() call.
    Ex:
    mock_db_session.output(
        [
            {'employee': 'John', 'Title': 'analyst'},
            {'employee': 'Mary', 'Title': 'PI'},
        ],
        repeat = True
    )


    If the test function performs two sqlalchemy queries, and we want to specify the return different entries
    for each query, (i.e., two queries returning 2 entries each), call mock_db_sssion.output() multiple times.
    Ex:
    mock_db_session.output(
        [
            {'employee': 'John', 'Title': 'analyst'},
            {'employee': 'Mary', 'Title': 'PI'},
        ]
    )
    mock_db_session.output(
        [
            {'employee': 'Bob', 'Title': 'scientist'},
            {'employee': 'Lisa', 'Title': 'researcher'},
        ]
    )
    Note that if repeat=True is specified in either mock_db_session.query call, the last query entry is repeated.

    Additionally, to intentionally raise a sqlalchemy.exc.SQLAlchemyError after a query is performed, set
    mock_db_session.raise_exception = True

    In order to check if the sqlalchemly add, commit, query and close was successfully called, we can checkthe
    mockup_session[key] boolean for True | False. The keys are:
    add, commit, query, close.

    To clear these keys, call mock_db_session.clear().
    """
    import sqlalchemy

    class MockQueryFields:
        def __init__(self, entries: dict):
            for key in entries:
                setattr(self, key, entries[key])

    class MockQueryResultIter:
        def __init__(self, mock_query: list):
            self.mock_query = mock_query
            self.idx = 0

        def __next__(self):
            retval = None
            if data_obj.limit_query:
                num_entries = data_obj.limit_query
            elif len(self.mock_query) > data_obj.entry_idx:
                num_entries = len(self.mock_query[data_obj.entry_idx])
            else:
                num_entries = 0

            if num_entries > self.idx:
                retval = self.mock_query[data_obj.entry_idx][self.idx]
            else:
                raise StopIteration

            self.idx += 1
            return retval

    class MockQueryResult:
        def __init__(self):
            self.entries = []

        def __iter__(self):
            return MockQueryResultIter(self.entries)

        def __getitem__(self, index: int):
            retval = None
            if len(self.entries) > data_obj.entry_idx and len(self.entries[data_obj.entry_idx]) > index:
                retval = self.entries[data_obj.entry_idx][index]
            return retval

        def add_mock_entry(self, entries: dict):
            sub_entries = []
            for entry in entries:
                sub_entries.append(MockQueryFields(entry))
                for key in entry:
                    setattr(self, key, entry[key])
            self.entries.append(sub_entries)

        def count(self, *args, **kwargs):
            val = 0
            if len(self.entries) > data_obj.entry_idx:
                num_entries = len(self.entries[data_obj.entry_idx])
                if data_obj.limit_query:
                    val = data_obj.limit_query if data_obj.limit_query <= num_entries else num_entries
                else:
                    val = num_entries
            return val

        def order_by(self, *args, **kwargs):
            return self

        def group_by(self, *args, **kwargs):
            return self

        def limit(self, number: int):
            data_obj.limit_query = number
            return self

        def all(self):
            return self

    class MockSessionQuery:
        result = MockQueryResult()

        @staticmethod
        def filter(*args, **kwargs):
            MockSessionQuery.result.entries = []
            data_obj.limit_query = 0
            data_obj.entry_idx += 1

            if data_obj.repeat_entry and data_obj.entry_idx >= len(data_obj.queries):
                data_obj.entry_idx = len(data_obj.queries) - 1

            for entry in data_obj.queries:
                MockSessionQuery.result.add_mock_entry(entry)
            return MockSessionQuery.result

        @staticmethod
        def filter_by(*args, **kwargs):
            return MockSessionQuery.filter(*args, **kwargs)

        @staticmethod
        def get(*args, **kwargs):
            return MockSessionQuery.filter(*args, **kwargs)

    class MockSession:
        @staticmethod
        def add(*args, **kwargs):
            data_obj.session['add'] = True
            if data_obj.raise_exception:
                raise sqlalchemy.exc.SQLAlchemyError()

        @staticmethod
        def commit(*args, **kwargs):
            data_obj.session['commit'] = True
            if data_obj.raise_exception:
                raise sqlalchemy.exc.SQLAlchemyError()

        @staticmethod
        def query(*args, **kwargs):
            data_obj.session['query'] = True
            if data_obj.raise_exception:
                raise sqlalchemy.exc.SQLAlchemyError()
            return MockSessionQuery

        @staticmethod
        def close(*args, **kwargs):
            data_obj.session['close'] = True

        @staticmethod
        def rollback(*args, **kwargs):
            data_obj.session['rollback'] = True

        @staticmethod
        def output(entries: list, repeat=False, raise_exception=False):
            data_obj.repeat_entry = repeat
            data_obj.raise_exception = raise_exception
            data_obj.queries.append(entries)

        @staticmethod
        def clear():
            data_obj.__init__()
            MockSessionQuery.result.entries = []

    @dataclass
    class Data:
        queries = []
        limit_query = 0
        entry_idx = -1
        repeat_entry = False
        raise_exception = False
        session = {
            'add': False,
            'commit': False,
            'query': False,
            'close': False,
            'rollback': False,
        }

    data_obj = Data()
    return MockSession()


@pytest.fixture
def mock_rpc_request(monkeypatch):
    """Fixture to mockup the rpc_client.RpcClient object for handling RMQ requests.

    To set the return value of the rpc_client.RpcClient.request() call, set the
    mockup_rpc_request.json = dictionary

    To intentionally throw an exception when calling the request, set the
    mockup_rpc_request.Exception = True and catch the rpc_client.ConnectionError exception.
    """
    from jaws_rpc import rpc_client

    class MockRpcClient():
        def __init__(self):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args, **kwargs):
            pass

        @staticmethod
        def request(*args, **kwargs):
            return data_obj.json

        @staticmethod
        def output(jsondata):
            data_obj.json = jsondata

    def mock_rpc_client(*args, **kwargs):
        if data_obj.Exception:
            raise rpc_client.ConnectionError("RPC client failed")
        return MockRpcClient()

    monkeypatch.setattr(rpc_client, 'RpcClient', mock_rpc_client)

    class Data():
        def __init__(self):
            self.json = None
            self.Exception = False

    data_obj = Data()
    return MockRpcClient()
