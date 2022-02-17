import pytest
import globus_sdk
from datetime import datetime

import jaws_central.analysis
import jaws_central.models_fsa
import jaws_central.config
from jaws_central.datatransfer_plugins import globus_transfer


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
        self.tag = "example"
        self.wdl_file = "example.wdl"
        self.json_file = "example.json"


@pytest.fixture()
def mock_database(monkeypatch):
    monkeypatch.setattr(jaws_central.models_fsa.db, "session", MockSession)


@pytest.fixture()
def mock_globus(monkeypatch):
    monkeypatch.setattr(globus_sdk, "NativeAppAuthClient", MockNativeAppAuthClient)
    monkeypatch.setattr(
        globus_sdk, "ConfidentialAppAuthClient", MockConfidentialAppAuthClient
    )
    monkeypatch.setattr(
        globus_sdk, "RefreshTokenAuthorizer", MockRefreshTokenAuthorizer
    )
    monkeypatch.setattr(
        globus_sdk, "ClientCredentialsAuthorizer", MockClientCredentialsAuthorizer
    )
    monkeypatch.setattr(globus_sdk, "TransferClient", MockGlobusTransferClient)


def test_list_sites(configuration):
    user = "test_user"
    result, code = jaws_central.analysis.list_sites(user)
    expected_sites = ["JGI", "NERSC"]
    assert isinstance(result, list)
    for site in result:
        assert isinstance(site, dict)
        assert "site_id" in site
        assert site["site_id"] in expected_sites
        assert "max_ram_gb" in site


def test_cancel_transfer(configuration, mock_database, mock_globus, mock_data_transfer):
    transfer_id = "without_error"
    status = jaws_central.analysis._cancel_transfer(mock_data_transfer, transfer_id)
    assert status == f"cancelled {transfer_id}"


def test_cancel_run(monkeypatch):
    """Test the cancel run functioning."""
    def get_run_cancelled(user_id, run_id):
        run = MockGetInactiveUploadRun()
        run.status = "cancelled"
        return run

    def get_run_regular(user_id, run_id):
        run = MockGetInactiveUploadRun()
        run.status = "running"
        return run

    def mock_cancel_run(user, run):
        run.status = "cancelled"
        run.result = "cancelled"
        return run

    def mock_cancel_transfer(transfer_task_id):
        pass

    def mock_rpc_call_cancel(user_id, run_id, method, params={}):
        assert isinstance(user_id, str)
        assert isinstance(run_id, int)
        assert method == "cancel_run"

    """Check if an exception is raised in case run is already cancelled"""
    monkeypatch.setattr(jaws_central.analysis, "rpc_call", mock_rpc_call_cancel)
    monkeypatch.setattr(jaws_central.analysis, "_get_run", get_run_cancelled)
    monkeypatch.setattr(jaws_central.analysis, "_cancel_run", mock_cancel_run)
    monkeypatch.setattr(jaws_central.analysis, "_cancel_transfer", mock_cancel_transfer)
    with pytest.raises(Exception):
        jaws_central.analysis.cancel_run("user", 123)

    """Check if no exception is raised in case run has a regular status"""
    monkeypatch.setattr(jaws_central.analysis, "_get_run", get_run_regular)
    jaws_central.analysis.cancel_run("user", 123)


def test_run_status_inactive_upload(monkeypatch):
    """ Tests returning run status with comments for a failed download status """

    def get_failed_run(user_id, run_id):
        return MockGetInactiveUploadRun()

    def mock_is_admin(user_id):
        return True

    monkeypatch.setattr(jaws_central.analysis, "_get_run", get_failed_run)
    monkeypatch.setattr(jaws_central.analysis, "_is_admin", mock_is_admin)

    comments = r"Globus transfer stalled."

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
    monkeypatch.setattr(jaws_central.analysis, "rpc_call", mock_rpc_call)
    jaws_central.analysis.run_metadata("test_user", 123)


def test_globus_transfer_path(configuration):
    globus = globus_transfer.DataTransfer()
    test_data = [
        (
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
            "/",
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
        ),  # noqa: E501
        (
            "/global/scratch/jaws/jaws-dev/uploads/mmelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
            "/global/scratch/jaws",
            "/jaws-dev/uploads/mmelara/7d2454c6-7052-4835-bb8d-e701c2d3df3e.wdl",
        ),  # noqa: E501
        (
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",  # noqa
            "/",
            "/global/cscratch1/sd/jaws/jaws-dev/users/mamelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",  # noqa
        ),  # noqa: E501
        (
            "/global/scratch/jaws/jaws-dev/uploads/mmelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",  # noqa
            "/global/scratch/jaws",
            "/jaws-dev/uploads/mmelara/CORI/global/u2/m/mamelara/jaws-quickstart-example/data/sample.fastq.bz2",
        ),  # noqa: E501
        (
            "/global/cfs/cdirs/jaws/data-repository-dev/mmelara/CORI/5247",
            "/",
            "/global/cfs/cdirs/jaws/data-repository-dev/mmelara/CORI/5247",
        ),  # noqa: E501
        (
            "/global/scratch/jaws/data-repository-dev/mmelara/JGI/3593",
            "/global/scratch/jaws",
            "/data-repository-dev/mmelara/JGI/3593",
        ),  # noqa: E501
    ]
    for (full_path, host_path, expected_path) in test_data:
        result = globus.virtual_transfer_path(full_path, host_path)
        assert result == expected_path


def test_run_info():

    test_run = jaws_central.models_fsa.Run(
        id=999,
        user_id="testuser",
        site_id="testsite",
        submission_id="testsubmissionid",
        input_site_id="testinputsite",
        input_endpoint="testinep",
        output_endpoint="testoutep",
        output_dir="testoutdir",
        wdl_file="./test.wdl",
        json_file="./test.json",
        tag="testtag",
        status="created",
        submitted=datetime.now(),
        updated=datetime.now(),
        cromwell_run_id=None,
        result=None,
        upload_task_id="testuptask",
        download_task_id="testdowntask",
    )
    expected_fields_complete = [
        "id",
        "submission_id",
        "cromwell_run_id",
        "result",
        "status",
        "status_detail",
        "site_id",
        "submitted",
        "updated",
        "input_site_id",
        "input_endpoint",
        "upload_task_id",
        "output_endpoint",
        "output_dir",
        "download_task_id",
        "user_id",
        "tag",
        "wdl_file",
        "json_file",
    ]
    expected_fields_partial = [
        "id",
        "result",
        "status",
        "status_detail",
        "site_id",
        "submitted",
        "updated",
        "input_site_id",
        "tag",
        "wdl_file",
        "json_file",
    ]

    complete_results = jaws_central.analysis._run_info(test_run, False, True)
    for key in complete_results:
        assert key in expected_fields_complete
    for key in expected_fields_complete:
        assert key in complete_results

    complete_results = jaws_central.analysis._run_info(test_run, True, False)
    for key in complete_results:
        assert key in expected_fields_complete
    for key in expected_fields_complete:
        assert key in complete_results

    partial_results = jaws_central.analysis._run_info(test_run, False, False)
    for key in partial_results:
        assert key in expected_fields_partial
    for key in expected_fields_partial:
        assert key in partial_results
