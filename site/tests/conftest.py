"""
File contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
import os
import shutil
import json
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass
from jaws_site import models, config, cromwell
import sqlalchemy
from sqlalchemy.orm.exc import NoResultFound
from moto import mock_s3
import boto3


S3_BUCKET = "site"
this_date = datetime.today()


@pytest.fixture()
def configuration(config_file):
    if config.conf is not None:
        config.Configuration._destructor()
    return config.Configuration(config_file)


@pytest.fixture()
def file_not_found_config(tmp_path):
    cfg = tmp_path / "wrong.ini"
    content = ""
    cfg.write_text(content)
    import os

    os.remove(cfg.as_posix())
    return cfg.as_posix()


@pytest.fixture
def config_file_zero(tmp_path):
    cfg = tmp_path / "jaws-site-empty.ini"
    content = ""
    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture
def config_file_wrong(tmp_path):
    cfg = tmp_path / "jaws-site-wrong.ini"
    content = """[RMQ]
    """
    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
    content = """
[RMQ]
host = localhost
vhost = jaws_test
user = jaws
password = password
num_threads = 5
max_retries = 3

[PERFORMANCE_METRICS]
done_dir = /tmp/done_dir
processed_dir = /tmp/processed_dir
running_dir = /tmp/running_dir
cleanup_time = 10

[GLOBUS]
client_id = AAAA
client_secret = BBBB
endpoint_id = jaws-testing
host_path = /global/scratch/jaws

[DB]
host = localhost
port = 3306
user = elmer_fudd
password = hunting
db = hunting_sites
dialect = mysql+mysqlconnector
host2 = ${JAWS_DB_HOST}
password2 = ${JAWS_DB_PASSWORD}123
test =

[SITE]
id = eagle
deployment = test
inputs_dir = /global/scratch/jaws/jaws-dev/inputs
max_user_active_runs = 1
max_transfer_threads = 10
file_permissions = 777

[CROMWELL]
url = http://localhost:8000

[AWS]
aws_access_key_id = AAAA
aws_secret_access_key = BBBB
s3_bucket = CCCC
"""
    cfg.write_text(content)
    yield cfg.as_posix()


@pytest.fixture
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
    content = """[RMQ]
host = https://rmq.nersc.gov
vhost = jaws_test
user = bugs_bunny
password = xqweasdasa
[PERFORMANCE_METRICS]
done_dir =
processed_dir =
running_dir =
cleanup_time =
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
deployment = prod
inputs_dir = /global/scratch/jaws/jaws-dev/inputs
[AWS]
aws_access_key_id = AAAA
aws_secret_access_key = BBBB
s3_bucket = CCCC
    """
    cfg.write_text(content)
    return cfg.as_posix()


def initRunModel(**kwargs):
    return models.Run(
        id=kwargs.get("id", "99"),
        user_id=kwargs.get("user_id", "test_user"),
        submission_id=kwargs.get("submission_id", "XXXX"),
        caching=(False if "caching" in kwargs and kwargs["caching"] is False else True),
        input_site_id=kwargs.get("input_site_id", "NERSC"),
        cromwell_run_id=kwargs.get("cromwell_run_id", "myid"),
        result=kwargs.get("result", "succeeded"),
        status=kwargs.get("status", "running"),
        submitted=kwargs.get("submitted", datetime.utcnow()),
        updated=kwargs.get("updated", datetime.utcnow()),
        workflow_root=kwargs.get("workflow_root", None),
    )


class MockRunModel:
    """Mock Run sqlalchemy orm model object with useable defaults."""

    def __init__(self, **kwargs):
        self.user_id = kwargs.get("user_id", "jaws")
        self.upload_task_id = kwargs.get("upload_task_id", "1")
        self.submission_id = kwargs.get("submission_id", "XXXX")
        self.cromwell_run_id = kwargs.get("cromwell_run_id", None)
        self.status = kwargs.get("status", "running")
        self.id = kwargs.get("id", "99")
        self.result = kwargs.get("result", None)
        self.output_endpoint = kwargs.get(
            "output_endpoint", "EXAMPLE_OUTPUT_ENDPOINT_ID"
        )
        self.output_dir = kwargs.get("output_dir", ".")
        self.download_task_id = kwargs.get("download_task_id", "325")
        self.email = kwargs.get("email", "jaws@vog.gov")
        self.cromwell_workflow_dir = (
            "/global/scratch/jaws/dev/cromwell-executions/test_wdl/myid"
        )
        self.input_site_id = kwargs.get("input_site_id", "CORI")
        self.caching = kwargs.get("caching", True)
        self.submitted = kwargs.get("submitted", datetime.utcnow())
        self.updated = kwargs.get("updated", datetime.utcnow())
        self.workflow_root = kwargs.get("workflow_root", None)


class MockRun:
    def __init__(self):
        self.data = MockRunModel()

    def metadata(self):
        return MockCromwellMetadata("localhost/api/workflows/v1", "xxx-xxx")

    def outputs(self, relpath=None):
        return {"test": "success"}

    def outfiles(self, complete=False, relpath=None):
        return {"test": "success"}

    def workflow_root(self):
        return {"test": "success"}

    def output_manifest(self, complete=False):
        return {"test": "success"}

    def cancel(self):
        return {"test": "success"}

    def errors(self):
        return {"test": "success"}

    def running_tasks(self):
        return {"test": "success"}

    def task_log(self):
        return {"test": "success"}

    def task_summary(self):
        return {"test": "success"}

    def mark_to_cancel(self):
        return {"test": "success"}


class MockTaskLog:
    def __init__(self, session, cromwell_run_id, logger=None):
        self.session = session
        self.cromwell_run_id = cromwell_run_id
        self.logger = logger
        self.data = []

    def table(self):
        return []


class MockCromwell:
    def __init__(self, url="localhost"):
        self.url = url
        self.workflows_url = f"{url}/api/workflows/v1"

    def get_metadata(self, workflow_id, data=None, cache={}):
        return MockCromwellMetadata(self.workflows_url, workflow_id, data, cache)

    def status(self):
        return True


class MockCromwellException:
    def __init__(self, url="localhost"):
        self.url = url
        self.workflows_url = f"{url}/api/workflows/v1"

    def status(self):
        raise Exception


class MockCromwellMetadata:
    def __init__(self, workflows_url, workflow_id, data=None, cache={}):
        self.workflows_url = workflows_url
        self.workflow_id = workflow_id
        self.tasks = None
        self.data = data

    def workflow_root(self):
        return "/example/cromwell-outputs/wdlName/workflowRoot"

    def outputs(self, **kwargs):
        return {}


class MockRpcClient:
    def __init__(self, params=None, logger=None):
        pass

    def request(self, method, params={}):
        response = {"result": None}
        return response


def mock_rpc_client(run):
    return MockRpcClient()


class MockResponses:
    """Requests response class"""

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


class MockTransferModel:
    """Mock Transfer sqlalchemy orm model object with useable defaults."""

    def __init__(self, **kwargs):
        self.manifest_json = "[]"
        if "manifest" in kwargs:
            assert type(kwargs["manifest"]) == list
            self.manifest_json = json.dumps(kwargs.get("manifest", []))
        elif "manifest_json" in kwargs:
            assert type(kwargs["manifest_json"]) == str
            self.manifest_json = kwargs.get("manifest_json", "[]")
        self.id = kwargs.get("id", "12")
        self.status = kwargs.get("status", "queued")
        self.submitted = kwargs.get("submitted", datetime.utcnow())
        self.updated = kwargs.get("updated", datetime.utcnow())
        self.src_base_dir = kwargs.get("src_base_dir", "/inputs")
        self.dest_base_dir = kwargs.get("dest_base_dir", "/inputs")
        self.src_site_id = kwargs.get("src_site_id", "NERSC")
        self.dest_site_id = kwargs.get("dest_site_id", "JGI")
        self.globus_transfer_id = kwargs.get("globus_transfer_id", None)
        self.reason = kwargs.get("reason", None)


class MockTransfer:
    def __init__(self, session, data, raise_exception=False):
        self.session = session
        self.data = data
        self.raise_exception = raise_exception

    @classmethod
    def from_id(cls, session, id):
        assert id is not None
        data = MockTransferModel(id=id)
        return cls(session, data)

    def status(self):
        return self.data.status

    def reason(self):
        return self.data.reason

    def submit_transfer(self):
        if self.raise_exception:
            raise Exception
        pass

    def cancel(self):
        self.data.status = "canceled"
        return self.data.reason


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
        CROMWELL_ID = "15774623-0f76-49ef-828c-3aa0ccd024f5"
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
def inputs_files():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "XXXX")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    for f in ["XXXX.wdl", "XXXX.json"]:
        file_path = os.path.join(root_dir, f)
        with open(file_path, "w") as outfile:
            outfile.write(f"output for {f}")
    file_path = os.path.join(root_dir, "XXXX.zip")
    with open(file_path, "wb") as outfile:
        outfile.write("output for XXXX.zip".encode())

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def inputs_files_without_zip():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "WWWW")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    for f in ["WWWW.wdl", "WWWW.json"]:
        file_path = os.path.join(root_dir, f)
        with open(file_path, "w") as outfile:
            outfile.write(f"output for {f}")

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def inputs_files_missing_json():
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "YYYY")
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    file_path = os.path.join(root_dir, "YYYY.wdl")
    with open(file_path, "w") as outfile:
        outfile.write("output for YYYY.wdl")

    yield

    shutil.rmtree(root_dir)


@pytest.fixture()
def inputs_files_empty_wdl():
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
    inputs_dir = tmp_path / "cwd/inputs/jaws"

    inputs_dir.mkdir(parents=True)

    for f in ["XXXX.wdl", "XXXX.json", "XXXX.orig.json", "XXXX.zip"]:
        file_path = inputs_dir / f
        file_path.write_text(f"output for {f}")

    yield

    shutil.rmtree(inputs_dir)


class MockSession:
    def __init__(self):
        self.needs_to_be_closed = False
        return

    def begin(self):
        return

    def begin_nested(self):
        return

    def commit(self):
        self.needs_to_be_closed = True
        return

    def close(self):
        self.needs_to_be_closed = False
        return

    def close_all(self):
        return

    def add(self, data):
        return

    def query(self, orm):
        return None

    def get(self, model, pk):
        return None


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

            if data_obj.status["succeeded"]:
                return Status.succeeded
            elif data_obj.status["failed"]:
                return Status.failed
            elif data_obj.status["transferring"]:
                return Status.transferring
            elif data_obj.status["inactive"]:
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
        return "globus_transfer"

    monkeypatch.setattr(DataTransferFactory, "__new__", mock_data_transfer)
    monkeypatch.setattr(Run, "_get_data_transfer_type", mock_get_data_transfer_type)

    @dataclass
    class DataTransferExceptions:
        DataTransferError = False
        DataTransferAPIError = False
        DataTransferNetworkError = False

    class Data:
        raises = DataTransferExceptions()
        status = {
            "succeeded": False,
            "failed": False,
            "transferring": False,
            "inactive": False,
        }

    data_obj = Data()
    return data_obj


@pytest.fixture
def mock_rpc_request(monkeypatch):
    """Fixture to mockup the rpc_client.RpcClient object for handling RMQ requests.

    To set the return value of the rpc_client.RpcClient.request() call, set the
    mockup_rpc_request.json = dictionary

    To intentionally throw an exception when calling the request, set the
    mockup_rpc_request.Exception = True and catch the rpc_client.ConnectionError exception.
    """
    from jaws_rpc import rpc_client

    class MockRpcClient:
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

    monkeypatch.setattr(rpc_client, "RpcClient", mock_rpc_client)

    class Data:
        def __init__(self):
            self.json = None
            self.Exception = False

    data_obj = Data()
    return MockRpcClient()


@pytest.fixture
def mock_sqlalchemy_session():
    """Fixture to mockup the sqlalchemy session.

    session = database.session()  # here, session is a sqlalchemy sessionmaker object.
    result = session.query(table).\
        filter_by(**filters).\
        order_by(table.start_date.desc()).\
        all()

    will be mocked up. The result variable can be set in the pytest function using the
    mock_sqlalchemy_session.query
    variable. This variable accepts a list of either dictionaries, or another list of dictionaries.

    Ex: if the test function performs one sqlalchemy query, and we want to specify the return value
    of that query, (i.e., one query returning two entries):
    mock_sqlalchemy_session.output(
        [
            {'employee': 'John', 'Title': 'analyst'},
            {'employee': 'Mary', 'Title': 'PI'},
        ]
    )

    If the test function performs multiple sqlalchemy queries, and we want to return the same result,
    add repeat=True to the mock_sqlalchemy_session.output() call.
    Ex:
    mock_sqlalchemy_session.output(
        [
            {'employee': 'John', 'Title': 'analyst'},
            {'employee': 'Mary', 'Title': 'PI'},
        ],
        repeat = True
    )


    If the test function performs two sqlalchemy queries, and we want to specify the return different
    entries for each query, (i.e., two queries returning 2 entries each), call
    mock_sqlalchemy_session.output() multiple times.
    Ex:
    mock_sqlalchemy_session.output(
        [
            {'employee': 'John', 'Title': 'analyst'},
            {'employee': 'Mary', 'Title': 'PI'},
        ]
    )
    mock_sqlalchemy_session.output(
        [
            {'employee': 'Bob', 'Title': 'scientist'},
            {'employee': 'Lisa', 'Title': 'researcher'},
        ]
    )
    Note that if repeat=True is specified in either mock_sqlalchemy_session.query call, the last
    query entry is repeated.

    Additionally, to intentionally raise a sqlalchemy.exc.SQLAlchemyError after a query is performed,
    set mock_sqlalchemy_session.raise_exception = True

    In order to check if the SQLAlchemy add, commit, query and close was successfully called,
    we can check the mockup_session[key] boolean for True | False. The keys are:
    add, commit, query, close.

    To clear these keys, call mock_sqlalchemy_session.clear().
    """

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
            if (
                len(self.entries) > data_obj.entry_idx
                and len(self.entries[data_obj.entry_idx]) > index
            ):
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
                    val = (
                        data_obj.limit_query
                        if data_obj.limit_query <= num_entries
                        else num_entries
                    )
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
            return self.entries[0] if len(self.entries) else []

        def one(self):
            if len(self.entries) > 0 and len(self.entries[0]) > 0:
                return self.entries[0][0]
            return None

        def one_or_none(self):
            if len(self.entries) > 0 and len(self.entries[0]) > 0:
                return self.entries[0][0]
            return None

        def filter(self, *args, **kwargs):
            return self

    class MockSessionQuery:
        result = MockQueryResult()

        @staticmethod
        def filter(*args, **kwargs):
            if data_obj.raise_exception_sqlalchemyerror:
                raise sqlalchemy.exc.SQLAlchemyError
            elif data_obj.raise_exception:
                raise Exception
            elif data_obj.raise_exception_integrityerror:
                raise sqlalchemy.exc.IntegrityError
            elif data_obj.raise_exception_noresultfound:
                raise NoResultFound

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
            if data_obj.raise_exception_sqlalchemyerror:
                raise sqlalchemy.exc.SQLAlchemyError
            elif data_obj.raise_exception:
                raise Exception
            elif data_obj.raise_exception_integrityerror:
                raise sqlalchemy.exc.IntegrityError

            return MockSessionQuery.filter(*args, **kwargs)

        @staticmethod
        def all(*args, **kwargs):
            return MockSessionQuery.filter(*args, **kwargs)

    class MockSessionNested:
        def __init__(self):
            self.data = data_obj

        def rollback(self):
            data_obj.session["rollback"] = True

    class MockSession:
        def __init__(self):
            self.data = data_obj

        def begin(self):
            return

        def begin_nested(self):
            return MockSessionNested()

        @staticmethod
        def add(*args, **kwargs):
            data_obj.session["add"] = True
            if data_obj.raise_exception_sqlalchemyerror:
                raise sqlalchemy.exc.SQLAlchemyError
            elif data_obj.raise_exception:
                raise Exception
            elif data_obj.raise_exception_integrityerror:
                raise sqlalchemy.exc.IntegrityError

        @staticmethod
        def commit(*args, **kwargs):
            data_obj.session["commit"] = True
            if data_obj.raise_exception_sqlalchemyerror:
                raise sqlalchemy.exc.SQLAlchemyError
            elif data_obj.raise_exception:
                raise Exception
            elif data_obj.raise_exception_integrityerror:
                raise sqlalchemy.exc.IntegrityError
            elif data_obj.raise_exception_commit:
                raise Exception

        @staticmethod
        def query(*args, **kwargs):
            data_obj.session["query"] = True
            if data_obj.raise_exception_integrityerror:
                raise sqlalchemy.exc.IntegrityError
            elif data_obj.raise_exception_sqlalchemyerror:
                raise sqlalchemy.exc.SQLAlchemyError
            elif data_obj.raise_exception:
                raise Exception

            return MockSessionQuery

        @staticmethod
        def close(*args, **kwargs):
            data_obj.session["close"] = True

        @staticmethod
        def rollback(*args, **kwargs):
            data_obj.session["rollback"] = True

        @staticmethod
        def output(entries: list, repeat=False, raise_exception=False):
            data_obj.repeat_entry = repeat
            data_obj.raise_exception = raise_exception
            data_obj.queries.append(entries)

        def clear(self):
            global data_obj
            data_obj = Data()
            data_obj.queries = []
            self.data = data_obj

    @dataclass
    class Data:
        queries = []
        limit_query = 0
        entry_idx = -1
        repeat_entry = False
        raise_exception = False
        raise_exception_commit = False
        raise_exception_sqlalchemyerror = False
        raise_exception_integrityerror = False
        raise_exception_noresultfound = False
        session = {
            "add": False,
            "commit": False,
            "query": False,
            "close": False,
            "rollback": False,
        }

    global data_obj
    data_obj = Data()
    return MockSession()


def initTransferModel(**kwargs):
    return models.Transfer(
        id=kwargs.get("id", "99"),
        status=kwargs.get("status", "created"),
        # src_site_id=kwargs.get("src_site_id", "NERSC"),
        src_base_dir=kwargs.get("src_base_dir", "/jaws-test/inputs"),
        # dest_site_id=kwargs.get("dest_site_id", "JGI"),
        dest_base_dir=kwargs.get("dest_base_dir", "/jaws-test/inputs"),
        manifest_json=kwargs.get("manifest_json", "{}"),
        # globus_transfer_id=kwargs.get("globus_transfer_id", None),
        reason=kwargs.get("reason", "reason"),
    )


@pytest.fixture
def mock_metadata(monkeypatch):
    class MockMetadata:
        def __init__(self):
            self.data = {
                "MOCK_METADATA": True,
                "workflowName": "unknown",
                "workflowRoot": "/data/cromwell-executions/example/ABCD",
                "status": "Running",
            }

        def started_running(self):
            return True

        def task_summary(self, last_attempts=False):
            return [
                {
                    "name": "test",
                    "shard_index": "shard_index",
                    "attempt": "attempt",
                    "cached": "cached",
                    "job_id": "123",
                    "execution_status": "succeeded",
                    "result": "result",
                    "failure_message": "failure_message",
                    "queue_start": "01-01-2022",
                    "run_start": "01-01-2022",
                    "run_end": "01-01-2022",
                    "queue_duration": "queue_duration",
                    "run_duration": "run_duration",
                    "call_root": "call_root",
                    "requested_time": "01-01-2022",
                    "requested_cpu": 15,
                    "requested_memory": "100",
                },
            ]

        def task_log(self):
            return {
                "name": "test",
                "cached": False,
                "execution_status": "succeeded",
                "queue_start": "01-01-2022",
                "run_start": "01-01-2022",
                "run_end": "01-01-2022",
                "queue_duration": "01-01-2022",
                "run_duration": "01-01-2022",
            }

        def get(self, param, default=None):
            return self.data.get(param, default)

        def workflow_root(self, executions_dir=None):
            return "/root"

        def errors(self):
            return "error"

        def running(self):
            return "running"

    def mock_cromwell_get_metadata(self, cromwell_run_id):
        assert type(cromwell_run_id) is str
        mock_metadata = MockMetadata()
        mock_metadata.data["runId"] = cromwell_run_id
        return mock_metadata

    monkeypatch.setattr(cromwell.Cromwell, "get_metadata", mock_cromwell_get_metadata)


@pytest.fixture
def s3():
    """Pytest fixture that creates the recipes bucket in
    the fake moto AWS account

    Yields a fake boto3 s3 client
    """
    with mock_s3():
        s3_client = boto3.client("s3")
        s3_client.create_bucket(
            Bucket=S3_BUCKET,
            CreateBucketConfiguration={"LocationConstraint": "us-west-1"},
        )
        yield s3_client


class MockLogger:
    def __init__(self):
        pass

    def error(self, message):
        pass

    def warning(self, message):
        pass

    def info(self, message):
        pass

    def debug(self, message):
        pass


@pytest.fixture
def setup_files(tmpdir):
    src_dir = tmpdir.mkdir("src")
    dst_dir = tmpdir.mkdir("dst")

    for i in range(1000):
        file = src_dir.join(f"file{i}.txt")
        file.write("content")

    yield str(src_dir), str(dst_dir)
