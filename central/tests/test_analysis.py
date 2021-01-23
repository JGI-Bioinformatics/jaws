import pytest
import globus_sdk
from datetime import datetime

import jaws_central.analysis
import jaws_central.models_fsa
import jaws_central.config


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-central.ini"
    content = """[JAWS]
name = jaws-dev
version = 2.0.1
docs_url = https://jaws-docs.readthedocs.io/en/latest/

[RPC_SERVER]
user = jaws
password = pppaass4
vhost = jaws_test

[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3306
user = jaws
password = passw0rd1
db = jaws

[GLOBUS]
client_id = AAAA
client_secret = BBBB

[SITE:LBNL]
host = rmq.jaws.gov
user = jaws
password = passw0rd2
vhost = jaws
queue = lbnl_rpc
globus_endpoint = XXXX
globus_basepath = "/global/scratch/jaws"
uploads_subdir = "uploads"
max_ram_gb = 1024

[SITE:NERSC]
host = rmq.jaws.gov
user = jaws
password = passw0rd2
vhost = jaws
queue = nersc_rpc
globus_endpoint = YYYY
globus_basepath = "/"
uploads_subdir = "/global/scratch/jaws/uploads"
max_ram_gb = 2048
"""
    cfg.write_text(content)
    return cfg.as_posix()


class MockClient:
    def __init__(self):
        pass


class MockAuthorizer:
    def __init__(self):
        pass


class MockNativeAppAuthClient:
    def __init__(self, client_id):
        pass


class MockConfidentialAppAuthClient:
    def __init__(self, client_id, client_secret):
        pass


class MockRefreshTokenAuthorizer:
    def __init__(self, token, client):
        pass


class MockClientCredentialsAuthorizer:
    def __init__(self, client, scopes):
        pass


class MockGlobusTransferClient:
    def __init__(self, authorizer=None):
        self.authorizer = authorizer

    def cancel_task(self, task_id):
        if task_id == "error":
            raise globus_sdk.GlobusAPIError("error")


class MockUser:
    @property
    def transfer_refresh_token(self):
        return "abcdefghijklmnopqrstuvwxyz"


class MockQuery:
    @staticmethod
    def get(user):
        return MockUser()


class MockSession:
    @staticmethod
    def query(user):
        return MockQuery()


class MockDb:
    @property
    def session(self):
        return MockSession()


class MockResponse:
    def __init__(self, status_code):
        self.status_code = status_code


class MockGetInactiveUploadRun:
    def __init__(self, *args, **kwargs):
        self.cromwell_run_id = "90b333ed-8095-474e-9ced-df86f3c99241"
        self.download_task_id = "5445cffc-963a-11ea-8ec9-02c81b96a709"
        self.id = 123
        self.input_endpoint = "9d6d994a-6d04-11e5-ba46-22000b92c6ec"
        self.input_site_id = "NERSC"
        self.output_dir = "/mydir/jaws-test-outdir"
        self.output_endpoint = "9d6d994a-6d04-11e5-ba46-22000b92c6ec"
        self.site_id = "NERSC"
        self.status = "upload inactive"
        self.result = None
        self.submission_id = "2ba6eb76-22e0-4d49-8eeb-b0c6683dfa30"
        self.submitted = datetime.strptime("2020-05-14 23:08:50", "%Y-%m-%d %H:%M:%S")
        self.updated = datetime.strptime("2020-05-14 23:27:15", "%Y-%m-%d %H:%M:%S")
        self.upload_task_id = "dfbdfb7a-9637-11ea-bf90-0e6cccbb0103"
        self.user_id = "dduck"


@pytest.fixture()
def mock_database(monkeypatch):
    monkeypatch.setattr(jaws_central.models_fsa.db, "session", MockSession)


@pytest.fixture()
def mock_globus(monkeypatch):
    monkeypatch.setattr(globus_sdk, "NativeAppAuthClient", MockNativeAppAuthClient)
    monkeypatch.setattr(globus_sdk, "ConfidentialAppAuthClient", MockConfidentialAppAuthClient)
    monkeypatch.setattr(
        globus_sdk, "RefreshTokenAuthorizer", MockRefreshTokenAuthorizer
    )
    monkeypatch.setattr(
        globus_sdk, "ClientCredentialsAuthorizer", MockClientCredentialsAuthorizer
    )
    monkeypatch.setattr(globus_sdk, "TransferClient", MockGlobusTransferClient)


@pytest.fixture()
def configuration(config_file):
    return jaws_central.config.Configuration(config_file)


def test_list_sites(configuration):
    user = "test_user"
    result, code = jaws_central.analysis.list_sites(user)
    assert isinstance(result, list)
    for site in result:
        assert isinstance(site, dict)
        assert 'site_id' in site
        assert 'max_ram_gb' in site


def test_cancel_transfer(configuration, mock_database, mock_globus):
    user = "test_user"
    transfer_id = "without_error"
    run_id = 99
    jaws_central.analysis._cancel_transfer(user, transfer_id, run_id)


def test_run_status_inactive_upload(monkeypatch):
    """ Tests returning run status with comments for a failed download status """

    def get_failed_run(user_id, run_id):
        return MockGetInactiveUploadRun()

    monkeypatch.setattr(jaws_central.analysis, "_get_run", get_failed_run)

    comments = r"Globus transfer stalled; try reactivating the endpoint. Please go to https://app.globus.org/file-manager, on the left side of the page, select ENDPOINTS, click the > icon to the right of the NERSC DTN endpoint, then click Activate."  # noqa"

    out_json, status = jaws_central.analysis.run_status("user", 123)
    assert out_json["status_detail"] == comments


def test_run_metadata(monkeypatch):

    def mock_get_run(user_id, run_id):
        return {}

    def mock_abort_if_pre_cromwell(run):
        return

    def mock_rpc_call(user_id, run_id, method, params={}):
        assert isinstance(user_id, str)
        assert isinstance(run_id, int)
        assert method == "run_metadata"

    monkeypatch.setattr(jaws_central.analysis, "_get_run", mock_get_run)
    monkeypatch.setattr(
        jaws_central.analysis, "_abort_if_pre_cromwell", mock_abort_if_pre_cromwell
    )
    monkeypatch.setattr(jaws_central.analysis, "_rpc_call", mock_rpc_call)
    jaws_central.analysis.run_metadata("test_user", 123)
