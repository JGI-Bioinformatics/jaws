import pytest
import os.path
import jaws_site
from jaws_site.runs import Run

# import jaws_site.globus
import tests.conftest
import globus_sdk


def test_run_constructor(monkeypatch):
    mock_session = tests.conftest.MockSession()
    run = Run(mock_session)
    assert run


def mock__update_run_status(self, new_status, timestamp):
    self.model.status = new_status
    self.model.updated = timestamp
    assert new_status is not None
    assert isinstance(new_status, str)
    assert len(new_status) > 0
    return


def mock__insert_run_log(self, status_from, status_to, timestamp, reason=None):
    assert isinstance(status_from, str)
    assert isinstance(status_to, str)


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


class MockCromwell:
    def __init__(self, url="localhost"):
        self.url = url
        self.workflows_url = f"{url}/api/workflows/v1"

    def get_metadata(self, workflow_id, data=None, cache={}):
        return MockCromwellMetadata(self.workflows_url, workflow_id, data, cache)


class MockCromwellMetadata:
    def __init__(self, workflows_url, workflow_id, data=None, cache={}):
        self.workflows_url = workflows_url
        self.workflow_id = workflow_id
        self.tasks = None
        self.data = data

    def workflow_root(self):
        return "/example/workflow/output"


class MockRpcClient:
    def __init__(self, params=None, logger=None):
        pass

    def request(self, method, params={}):
        response = {"result": None}
        return response


def mock_rpc_client(run):
    return MockRpcClient()


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
    mock_session = tests.conftest.MockSession()
    run = Run(mock_session)
    proc = run.operations.get(status, None)
    assert callable(proc)


@pytest.mark.parametrize(
    "statuses",
    [{"status": "FAILED"}, {"status": "INACTIVE"}],
    ids=["failed", "inactive"],
)
def test_check_if_upload_complete(statuses, monkeypatch):
    """
    Tests check_if_upload_complete from Run class. This only tests two
    statuses since the 'SUCCEEDED' status calls another method we want to
    test later and separately.
    """
    mock_session = tests.conftest.MockSession()
    model = tests.conftest.MockRunModel(status="uploading")
    run = Run(mock_session, model=model)

    def mock_globus_service(run):
        return MockGlobusService(statuses)

    monkeypatch.setattr(jaws_site.globus, "GlobusService", mock_globus_service)
    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)
    run.check_if_upload_complete()


def test_file_size(uploads_files_empty_wdl):

    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "ZZZZ")

    test_files = [
        (os.path.join(root_dir, "ZZZZ.json"), 2),
        (os.path.join(root_dir, "ZZZZ.wdl"), 0),
    ]

    for (test_file, expected) in test_files:
        test_size = jaws_site.runs.file_size(test_file)
        assert test_size == expected

    test_size = jaws_site.runs.file_size(os.path.join(root_dir, "ZZZZ.zip"))
    assert test_size is None


def test_get_run_input(
    monkeypatch, uploads_files, uploads_files_missing_json, uploads_files_without_zip
):
    def mock_uploads_file_path(self):
        submission_id = self.model.submission_id
        home_dir = os.path.expanduser("~")
        root_dir = os.path.join(home_dir, submission_id)
        return os.path.join(root_dir, submission_id)

    monkeypatch.setattr(jaws_site.runs.Run, "uploads_file_path", mock_uploads_file_path)

    mock_session = tests.conftest.MockSession()

    # test 1: valid input
    model1 = tests.conftest.MockRunModel(status="upload complete", submission_id="XXXX")
    run1 = Run(mock_session, model=model1)
    infiles = run1.get_run_input()
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "XXXX")
    assert infiles[0] == os.path.join(root_dir, "XXXX.wdl")
    assert infiles[1] == os.path.join(root_dir, "XXXX.json")
    assert infiles[2] == os.path.join(root_dir, "XXXX.zip")
    assert infiles[3] is None

    # test 2: invalid input
    model2 = tests.conftest.MockRunModel(status="upload complete", submission_id="YYYY")
    run2 = Run(mock_session, model=model2)
    with pytest.raises(jaws_site.runs.DataError):
        infiles = run2.get_run_input()

    # test 3: valid input, no subworkflows zip
    model3 = tests.conftest.MockRunModel(status="upload complete", submission_id="WWWW")
    run3 = Run(mock_session, model=model3)
    infiles = run3.get_run_input()
    home_dir = os.path.expanduser("~")
    root_dir = os.path.join(home_dir, "WWWW")
    assert infiles[0] == os.path.join(root_dir, "WWWW.wdl")
    assert infiles[1] == os.path.join(root_dir, "WWWW.json")
    assert infiles[2] is None
    assert infiles[3] is None


def test_submit_run(monkeypatch, uploads_files):
    def mock_cromwell_submit(self, wdl_file, json_file, zip_file):
        return "ABCD-EFGH"

    mock_session = tests.conftest.MockSession()
    mock_model = tests.conftest.MockRunModel(status="uploading")
    run = Run(mock_session, model=mock_model)

    monkeypatch.setattr(jaws_site.cromwell.Cromwell, "submit", mock_cromwell_submit)
    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)

    run.submit_run()


@pytest.fixture
def mock_path(tmp_path):
    def tmp_path_mock(*args, **kwargs):
        return (tmp_path / "cwd/uploads/jaws/XXXX").as_posix()

    return tmp_path_mock


def test_transfer_results(monkeypatch, transfer_dirs, mock_path, tmp_path):
    def mock_cp_infile_to_outdir(self, path, suffix, dest, required=True):
        pass

    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)
    monkeypatch.setattr(Run, "uploads_file_path", mock_path)
    monkeypatch.setattr(Run, "_cp_infile_to_outdir", mock_cp_infile_to_outdir)

    mock_session = tests.conftest.MockSession()
    mock_model = tests.conftest.MockRunModel(
        status="completed", cromwell_run_id="EXAMPLE_CROMWELL_ID"
    )

    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.runs, "globus", MockGlobusService({"status": "SUCCEEDED"})
    )
    monkeypatch.setattr(jaws_site.runs, "cromwell", MockCromwell())

    run.transfer_results()

    assert run.model.download_task_id == "325"


def test_failed_transfer_result(monkeypatch, transfer_dirs, mock_path, tmp_path):
    def mock_cp_infile_to_outdir(self, path, suffix, dest, required=True):
        pass

    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)
    monkeypatch.setattr(Run, "uploads_file_path", mock_path)
    monkeypatch.setattr(Run, "_cp_infile_to_outdir", mock_cp_infile_to_outdir)

    mock_session = tests.conftest.MockSession()
    mock_model = tests.conftest.MockRunModel(
        status="completed", cromwell_run_id="EXAMPLE_CROMWELL_ID"
    )

    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.runs, "globus", MockGlobusWithError({"status": "FAILED"})
    )
    monkeypatch.setattr(jaws_site.runs, "cromwell", MockCromwell())

    with pytest.raises(globus_sdk.GlobusError):
        run.transfer_results()


@pytest.mark.parametrize(
    "status,expected",
    [
        ({"status": "SUCCEEDED"}, "download complete"),
        ({"status": "FAILED"}, "download failed"),
    ],
)
def test_check_if_download_complete(status, expected, monkeypatch):

    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)

    mock_session = tests.conftest.MockSession()
    mock_model = tests.conftest.MockRunModel(status="downloading")
    run = Run(mock_session, model=mock_model)

    monkeypatch.setattr(jaws_site.runs, "globus", MockGlobusService(status))

    run.check_if_download_complete()

    assert run.status == expected


def test_check_run_cromwell_status(monkeypatch):

    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)

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

    mock_session = tests.conftest.MockSession()

    # test: submitted -> queued
    mock_model = tests.conftest.MockRunModel(status="queued")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.tasks, "get_run_status", mock_get_run_status_queued)
    run.check_run_cromwell_status()
    assert run.status == "queued"

    # test: queued -> queued
    mock_model = tests.conftest.MockRunModel(status="queued")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.tasks, "get_run_status", mock_get_run_status_none)
    run.check_run_cromwell_status()
    assert run.status == "queued"

    # test: queued -> running
    mock_model = tests.conftest.MockRunModel(status="queued")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.tasks, "get_run_status", mock_get_run_status_running)
    run.check_run_cromwell_status()
    assert run.status == "running"

    # test: queued -> succeeded
    mock_model = tests.conftest.MockRunModel(status="queued")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    run.check_run_cromwell_status()
    assert run.status == "succeeded"

    # test: queued -> failed
    mock_model = tests.conftest.MockRunModel(status="queued")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    run.check_run_cromwell_status()
    assert run.status == "failed"

    # test: running -> succeeded
    mock_model = tests.conftest.MockRunModel(status="running")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    run.check_run_cromwell_status()
    assert run.status == "succeeded"

    # test: running -> failed
    mock_model = tests.conftest.MockRunModel(status="running")
    run = Run(mock_session, model=mock_model)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    run.check_run_cromwell_status()
    assert run.status == "failed"
