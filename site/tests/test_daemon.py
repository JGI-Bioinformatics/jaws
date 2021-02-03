import pytest
import globus_sdk

from jaws_site.daemon import Daemon
import jaws_site.globus
import tests.conftest
import requests


def mock_update_run_status(self, run, new_status, reason=None):
    run.status = new_status
    assert isinstance(new_status, str)
    return


class MockGlobusService:
    def __init__(self, status, transfer_result={"task_id": "325"}):
        self.status = status
        self.transfer_result = transfer_result

    def transfer_status(self, task_id):
        return self.status['status']

    def submit_transfer(self, run_id, endpoint, src_dir, dest_dir):
        return self.transfer_result["task_id"]


class MockGlobusWithError(MockGlobusService):

    def submit_transfer(self, run_id, endpoint, src_dir, dest_dir):
        raise globus_sdk.GlobusError()


@pytest.mark.parametrize(
    "status",
    [
        "uploading",
        "upload complete",
        "submitted",
        "queued",
        "running",
        "failed",
        "succeeded",
        "downloading",
    ],
)
def test_check_operations_table(status):
    """Check one of the many possible entries in the operations table, which should return a method."""
    daemon = Daemon()
    proc = daemon.operations.get(status, None)
    assert callable(proc)


@pytest.mark.parametrize(
    "statuses",
    [{"status": "FAILED"}, {"status": "INACTIVE"}],
    ids=["failed", "inactive"],
)
def test_check_if_upload_complete(statuses, monkeypatch):
    """
    Tests check_if_upload_complete from Daemon class. This only tests two
    statuses since the 'SUCCEEDED' status calls another method we want to
    test later and separately.
    """
    daemon = Daemon()
    run = tests.conftest.MockRun(status="uploading")

    def mock_globus_service(daemon):
        return MockGlobusService(statuses)

    monkeypatch.setattr(jaws_site.globus, "GlobusService", mock_globus_service)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    daemon.check_if_upload_complete(run)


def test_submit_run(monkeypatch, uploads_files):
    def workflows_post(url, files={}):
        return tests.conftest.MockResponses({"id": "2"}, 201)

    daemon = Daemon()
    run = tests.conftest.MockRun(status="upload complete")

    monkeypatch.setattr(requests, "post", workflows_post)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)

    daemon.submit_run(run)


@pytest.fixture
def mock_path(tmp_path):
    def tmp_path_mock(*args, **kwargs):
        return (tmp_path / "cwd/uploads/jaws/2").as_posix()
    return tmp_path_mock


def test_transfer_results(monkeypatch, transfer_dirs, mock_path, tmp_path):

    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    monkeypatch.setattr(Daemon, "get_uploads_file_path", mock_path)

    daemon = Daemon()
    monkeypatch.setattr(daemon, "globus", MockGlobusService({"status": "SUCCEEDED"}))
    monkeypatch.setattr(daemon, "session", tests.conftest.MockSession())

    run = tests.conftest.MockRun(status="running", cromwell_run_id="EXAMPLE_CROMWELL_ID")
    run.cromwell_workflow_dir = tmp_path.as_posix()

    daemon.transfer_results(run)

    assert run.download_task_id
    assert run.download_task_id == "325"


def test_failed_transfer_result(monkeypatch, transfer_dirs, mock_path, tmp_path):
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    monkeypatch.setattr(Daemon, "get_uploads_file_path", mock_path)

    daemon = Daemon()
    monkeypatch.setattr(daemon, "globus", MockGlobusWithError({"status": "FAILED"}))
    monkeypatch.setattr(daemon, "session", tests.conftest.MockSession())

    run = tests.conftest.MockRun(status="running", cromwell_run_id="EXAMPLE_CROMWELL_ID")
    run.cromwell_workflow_dir = tmp_path.as_posix()

    with pytest.raises(globus_sdk.GlobusError):
        daemon.transfer_results(run)


@pytest.mark.parametrize(
    "status,expected",
    [({"status": "SUCCEEDED"}, "download complete"), ({"status": "FAILED"}, "download failed")],
)
def test_check_if_download_complete(status, expected, monkeypatch):

    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    daemon = Daemon()
    monkeypatch.setattr(daemon, "globus", MockGlobusService(status))

    run = tests.conftest.MockRun(status="downloading")
    daemon.check_if_download_complete(run)

    assert run.status == expected
