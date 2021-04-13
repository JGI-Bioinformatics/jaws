import pytest
import globus_sdk
import os.path

from jaws_site.daemon import Daemon
import jaws_site.globus
import tests.conftest


def mock_update_run_status(self, run, new_status, reason=None):
    run.status = new_status
    assert isinstance(new_status, str)
    return


class MockGlobusService:
    def __init__(self, status, transfer_result={"task_id": "325"}):
        self.status = status
        self.transfer_result = transfer_result

    def transfer_status(self, task_id):
        return self.status["status"]

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


def test_file_size(uploads_files_empty_wdl):

    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "ZZZZ")

    test_files = [
        (os.path.join(root_dir, "ZZZZ.json"), 2),
        (os.path.join(root_dir, "ZZZZ.wdl"), 0),
    ]

    for (test_file, expected) in test_files:
        test_size = jaws_site.daemon.file_size(test_file)
        assert test_size == expected

    test_size = jaws_site.daemon.file_size(os.path.join(root_dir, "ZZZZ.zip"))
    assert test_size is None


def test_get_run_input(
    monkeypatch, uploads_files, uploads_files_missing_json, uploads_files_without_zip
):
    def mock_get_uploads_file_path(self, run):
        submission_id = run.submission_id
        home_dir = os.path.expanduser("~")
        root_dir = os.path.join(home_dir, submission_id)
        return os.path.join(root_dir, submission_id)

    monkeypatch.setattr(
        jaws_site.daemon.Daemon, "get_uploads_file_path", mock_get_uploads_file_path
    )

    daemon = Daemon()

    # test 1: valid input
    run1 = tests.conftest.MockRun(status="upload complete", submission_id="XXXX")
    infiles = daemon.get_run_input(run1)
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "XXXX")
    assert infiles[0] == os.path.join(root_dir, "XXXX.wdl")
    assert infiles[1] == os.path.join(root_dir, "XXXX.json")
    assert infiles[2] == os.path.join(root_dir, "XXXX.zip")

    # test 2: invalid input
    run2 = tests.conftest.MockRun(status="upload complete", submission_id="YYYY")
    with pytest.raises(jaws_site.daemon.DataError):
        infiles = daemon.get_run_input(run2)

    # test 3: valid input, no subworkflows zip
    run1 = tests.conftest.MockRun(status="upload complete", submission_id="WWWW")
    infiles = daemon.get_run_input(run1)
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "WWWW")
    assert infiles[0] == os.path.join(root_dir, "WWWW.wdl")
    assert infiles[1] == os.path.join(root_dir, "WWWW.json")
    assert len(infiles) == 2


def test_submit_run(monkeypatch, uploads_files):
    def mock_cromwell_submit(self, wdl_file, json_file, zip_file):
        return "ABCD-EFGH"

    daemon = Daemon()
    run = tests.conftest.MockRun(status="upload complete")

    monkeypatch.setattr(jaws_site.cromwell.Cromwell, "submit", mock_cromwell_submit)
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)

    daemon.submit_run(run)


@pytest.fixture
def mock_path(tmp_path):
    def tmp_path_mock(*args, **kwargs):
        return (tmp_path / "cwd/uploads/jaws/XXXX").as_posix()

    return tmp_path_mock


def test_transfer_results(monkeypatch, transfer_dirs, mock_path, tmp_path):

    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    monkeypatch.setattr(Daemon, "get_uploads_file_path", mock_path)

    daemon = Daemon()
    monkeypatch.setattr(daemon, "globus", MockGlobusService({"status": "SUCCEEDED"}))
    monkeypatch.setattr(daemon, "session", tests.conftest.MockSession())

    run = tests.conftest.MockRun(
        status="running", cromwell_run_id="EXAMPLE_CROMWELL_ID"
    )
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

    run = tests.conftest.MockRun(
        status="running", cromwell_run_id="EXAMPLE_CROMWELL_ID"
    )
    run.cromwell_workflow_dir = tmp_path.as_posix()

    with pytest.raises(globus_sdk.GlobusError):
        daemon.transfer_results(run)


@pytest.mark.parametrize(
    "status,expected",
    [
        ({"status": "SUCCEEDED"}, "download complete"),
        ({"status": "FAILED"}, "download failed"),
    ],
)
def test_check_if_download_complete(status, expected, monkeypatch):

    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)
    daemon = Daemon()
    monkeypatch.setattr(daemon, "globus", MockGlobusService(status))

    run = tests.conftest.MockRun(status="downloading")
    daemon.check_if_download_complete(run)

    assert run.status == expected


def test_check_run_cromwell_status(monkeypatch):
    def mock_get_status_running(self, run_id):
        return "Running"

    def mock_get_status_succeeded(self, run_id):
        return "Succeeded"

    def mock_get_status_failed(self, run_id):
        return "Failed"

    def mock_get_run_status_queued(session, run_id):
        return "queued"

    def mock_get_run_status_none(session, run_id):
        return None

    def mock_get_run_status_running(session, run_id):
        return "running"

    daemon = Daemon()
    monkeypatch.setattr(Daemon, "update_run_status", mock_update_run_status)

    # test: submitted -> queued
    run = tests.conftest.MockRun(status="submitted")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.tasks, "get_run_status", mock_get_run_status_queued)
    daemon.check_run_cromwell_status(run)
    assert run.status == "queued"

    # test: queued -> queued
    run = tests.conftest.MockRun(status="queued")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.tasks, "get_run_status", mock_get_run_status_none)
    daemon.check_run_cromwell_status(run)
    assert run.status == "queued"

    # test: queued -> running
    run = tests.conftest.MockRun(status="queued")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.tasks, "get_run_status", mock_get_run_status_running)
    daemon.check_run_cromwell_status(run)
    assert run.status == "running"

    # test: queued -> succeeded
    run = tests.conftest.MockRun(status="queued")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    daemon.check_run_cromwell_status(run)
    assert run.status == "succeeded"

    # test: queued -> failed
    run = tests.conftest.MockRun(status="queued")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    daemon.check_run_cromwell_status(run)
    assert run.status == "failed"

    # test: running -> succeeded
    run = tests.conftest.MockRun(status="running")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    daemon.check_run_cromwell_status(run)
    assert run.status == "succeeded"

    # test: running -> failed
    run = tests.conftest.MockRun(status="running")
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    daemon.check_run_cromwell_status(run)
    assert run.status == "failed"
