import pytest
import globus_sdk
import pathlib
import os

from jaws_site.jawsd import Daemon
import tests.conftest
import requests


@pytest.mark.parametrize("statuses", [{"status": "FAILED"},
                                      {"status": "INACTIVE"}],
                         ids=["failed", "inactive"])
def test_check_transfer_status(statuses, monkeypatch):
    """
    Tests check_transfer_status from Daemon class. This only tests two
    statuses since the 'SUCCEEDED' status calls another method we want to
    test later and separately.
    """
    monkeypatch.setattr(Daemon, "_query_user_id", tests.conftest.query_jaws_id)
    db = tests.conftest.MockDb()
    jawsd = Daemon(db)
    run = tests.conftest.MockRun("jaws", "1", "2", "3", "myid", "submitted", "mydir")

    def mock_authorize_client(jawsd, token):
        return tests.conftest.MockTransferClient(statuses)

    monkeypatch.setattr(Daemon,
                        "_authorize_transfer_client",
                        mock_authorize_client)
    jawsd.check_transfer_status(run)


def test_submit_run(monkeypatch, staging_files):

    def workflows_post(url, files={}):
        return tests.conftest.MockResponses({"id": "2"}, 201)

    db = tests.conftest.MockDb()
    jawsd = Daemon(db)
    run = tests.conftest.MockRun("jaws", "1", "2", "3", "myid", "submitted", "mydir")

    monkeypatch.setattr(requests, "post", workflows_post)

    jawsd.submit_run(run)


test_data = [({"status": "Running"}, "submitted", "running"),
             ({"status": "Failed"}, "running", "failed"),
             ({"status": "Aborted"}, "aborting", "aborted")]


@pytest.mark.parametrize("status,current_status,expected_run", test_data)
def test_check_run_status(status, current_status, expected_run, monkeypatch):

    def workflows_get(url):
        return tests.conftest.MockResponses(status, 200)

    def prepare_failed_run_output(run):
        return

    db = tests.conftest.MockDb()
    jawsd = Daemon(db)
    run = tests.conftest.MockRun("jaws", "1", "2", "3", "myid", current_status, "mydir")

    monkeypatch.setattr(requests, "get", workflows_get)
    monkeypatch.setattr(jawsd, "prepare_failed_run_output", prepare_failed_run_output)

    jawsd.check_run_status(run)
    print(f"run.status={run.status}")
    print(f"expected_run={expected_run}")

    assert run.status == expected_run


def test_prepare_output():
    # TODO: Need to implement test for this
    pass


@pytest.mark.parametrize("status,current_status,expected_run", test_data)
def test_prepare_failed_run_output(status, current_status, expected_run, monkeypatch, cromwell_run_dir, user_dir):

    def workflows_get(url):
        return tests.conftest.MockResponses({"workflowRoot": cromwell_run_dir}, 200)

    def mock_authorize_client(jawd, token):
        return tests.conftest.MockTransferClientWithCopy({"status": "running"})

    def mock_user_query(jawsd, token):
        return tests.conftest.MockUserQuery(token)

    monkeypatch.setattr(requests, "get", workflows_get)
    monkeypatch.setattr(Daemon, "_authorize_transfer_client",
                        mock_authorize_client)
    monkeypatch.setattr(Daemon, "_query_user_id",
                        mock_user_query)
    monkeypatch.setattr(globus_sdk, "TransferData",
                        tests.conftest.MockTransferDataWithCopy)

    db = tests.conftest.MockDb()
    jawsd = Daemon(db)
    jawsd.results_dir = pathlib.Path(cromwell_run_dir).parent
    run = tests.conftest.MockRun("jaws", "1", "2", "3", "jaws_dir", "submitted", user_dir)
    jawsd.prepare_failed_run_output(run)

    for fname in os.listdir(cromwell_run_dir):
        bname = os.path.basename(fname)
        assert os.path.isdir(os.path.join(user_dir, bname))


def test_transfer_results(mock_query_user_id, monkeypatch):

    def mock_authorize_client(jawd, token):
        return tests.conftest.MockTransferClient({"status": "running"})

    monkeypatch.setattr(Daemon, "_authorize_transfer_client",
                        mock_authorize_client)
    monkeypatch.setattr(globus_sdk, "TransferData",
                        tests.conftest.MockTransferData)

    db = tests.conftest.MockDb()
    jawsd = Daemon(db)
    run = tests.conftest.MockRun("jaws", "1", "2", "3", "myid", "running", "mydir")

    jawsd.transfer_results(run)

    assert run.download_task_id
    assert run.download_task_id == '325'


@pytest.mark.parametrize("status,expected", [({"status": "SUCCEEDED"},
                                              "finished"),
                                             ({"status": "FAILED"},
                                              "download failed")])
def test_check_if_downloads_completes(status, expected, monkeypatch,
                                      mock_query_user_id):

    def mock_authorize_client(jawd, token):
        return tests.conftest.MockTransferClient(status)

    monkeypatch.setattr(Daemon, "_authorize_transfer_client",
                        mock_authorize_client)

    db = tests.conftest.MockDb()
    jawsd = Daemon(db)
    run = tests.conftest.MockRun("jaws", "1", "2", "3", "myid", "submitted", "mydir")

    jawsd.check_if_download_complete(run)

    assert run.status == expected
