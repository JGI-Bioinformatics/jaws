import pytest
import globus_sdk

from jaws_site.daemon import Daemon
import tests.conftest
import requests


def mock_update_run_status(self, run, new_status, reason=None):
    run.status = new_status
    assert isinstance(new_status, str)
    return


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
        "ready",
        "downloading",
    ],
)
def test_check_operations_table(status):
    """Check one of the many possible entries in the operations table, which should return a method."""
    jawsd = Daemon()
    proc = jawsd.operations.get(status, None)
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
    jawsd = Daemon()
    run = tests.conftest.MockRun(status="uploading")

    def mock_authorize_client(jawsd, token):
        return tests.conftest.MockTransferClient(statuses)

    monkeypatch.setattr(Daemon, "_authorize_transfer_client", mock_authorize_client)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    jawsd.check_if_upload_complete(run)


def test_submit_run(monkeypatch, staging_files):
    def workflows_post(url, files={}):
        return tests.conftest.MockResponses({"id": "2"}, 201)

    jawsd = Daemon()
    run = tests.conftest.MockRun(status="upload complete")

    monkeypatch.setattr(requests, "post", workflows_post)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)

    jawsd.submit_run(run)


def test_transfer_results(monkeypatch):
    def mock_authorize_client(jawd, token):
        return tests.conftest.MockTransferClient({"status": "running"})

    monkeypatch.setattr(Daemon, "_authorize_transfer_client", mock_authorize_client)
    monkeypatch.setattr(globus_sdk, "TransferData", tests.conftest.MockTransferData)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)

    jawsd = Daemon()
    run = tests.conftest.MockRun(status="running", cromwell_run_id="EXAMPLE_CROMWELL_ID")

    jawsd.transfer_results(run)

    assert run.download_task_id
    assert run.download_task_id == "325"


@pytest.mark.parametrize(
    "status,expected",
    [({"status": "SUCCEEDED"}, "download complete"), ({"status": "FAILED"}, "download failed")],
)
def test_check_if_download_complete(status, expected, monkeypatch):
    def mock_authorize_client(jawd, token):
        return tests.conftest.MockTransferClient(status)

    def mock_get_globus_transfer_status(self, run, task_id):
        return status["status"]

    monkeypatch.setattr(Daemon, "_authorize_transfer_client", mock_authorize_client)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    monkeypatch.setattr(
        Daemon, "_get_globus_transfer_status", mock_get_globus_transfer_status
    )

    jawsd = Daemon()
    run = tests.conftest.MockRun(status="downloading")

    jawsd.check_if_download_complete(run)

    assert run.status == expected
