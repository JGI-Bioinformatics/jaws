import pytest
import globus_sdk

import jaws_central.analysis
import jaws_central.models
import jaws_central.config


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


class MockClient:
    def __init__(self):
        pass


class MockAuthorizer:
    def __init__(self):
        pass


class MockNativeAppAuthClient:
    def __init__(self, client_id):
        pass


class MockRefreshTokenAuthorizer:
    def __init__(self, token, client):
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


@pytest.fixture()
def mock_database(monkeypatch):
    monkeypatch.setattr(jaws_central.models.db, "session", MockSession)


@pytest.fixture()
def mock_globus(monkeypatch):
    monkeypatch.setattr(globus_sdk, "NativeAppAuthClient", MockNativeAppAuthClient)
    monkeypatch.setattr(globus_sdk, "RefreshTokenAuthorizer", MockRefreshTokenAuthorizer)
    monkeypatch.setattr(globus_sdk, "TransferClient", MockGlobusTransferClient)


@pytest.fixture()
def configuration(config_file):
    return jaws_central.config.Configuration(config_file)


def test_cancel_transfer(configuration, mock_database, mock_globus):
    user = "test_user"
    transfer_id = "without_error"
    run_id = 99
    jaws_central.analysis._cancel_transfer(user, transfer_id, run_id)
