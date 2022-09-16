import pytest
import json
from deepdiff import DeepDiff
import jaws_site
from jaws_site.runs import (
    Run,
    RunLog,
    RunDbError,
    RunNotFoundError,
    RunFileNotFoundError,
)
from tests.conftest import MockSession, MockRunModel, this_date, initRunModel
from jaws_site.cromwell import Cromwell, CromwellError
import io


def mock__update_run_status(self, new_status, timestamp):
    self.data.status = new_status
    self.data.updated = timestamp
    assert new_status is not None
    assert isinstance(new_status, str)
    assert len(new_status) > 0
    return


def mock__insert_run_log(self, status_from, status_to, timestamp, reason=None):
    assert isinstance(status_from, str)
    assert isinstance(status_to, str)


def test_constructor():
    mock_session = MockSession()
    mock_data = MockRunModel()
    run = Run(mock_session, mock_data)
    assert run


def test_status():
    mock_session = MockSession()
    mock_data = MockRunModel(status="queued")
    run = Run(mock_session, mock_data)
    assert run.status() == "queued"


@pytest.mark.parametrize(
    "status",
    [
        "ready",
        "submitted",
        "queued",
        "running",
        "failed",
        "succeeded",
    ],
)
def test_check_operations_table(status):
    """Check one of the many possible entries in the operations table, which should return a method."""
    mock_session = MockSession()
    mock_data = MockRunModel()
    run = Run(mock_session, mock_data)
    proc = run.operations.get(status, None)
    assert callable(proc)


def test_mark_to_cancel(monkeypatch):
    def mock_update_run_status(self, new_status):
        assert new_status is not None
        assert type(new_status) == str
        assert new_status != self.data.status
        self.data.status = new_status

    monkeypatch.setattr(Run, "update_run_status", mock_update_run_status)

    mock_session = MockSession()

    # test1: run hasn't been cancelled yet
    mock_data = MockRunModel(status="queued", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.status() == "cancel"

    # test2: run was already marked to be cancelled
    mock_data = MockRunModel(status="cancel", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.status() == "cancel"

    # test3: run was already cancelled
    mock_data = MockRunModel(status="cancelled", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.status() == "cancelled"


def test_cancel(monkeypatch):
    def mock_cromwell_abort(self, cromwell_run_id):
        assert type(cromwell_run_id) is str

    def mock_update_run_status(self, new_status):
        assert new_status is not None
        assert type(new_status) == str
        assert new_status != self.data.status
        self.data.status = new_status

    monkeypatch.setattr(jaws_site.cromwell.Cromwell, "abort", mock_cromwell_abort)
    monkeypatch.setattr(Run, "update_run_status", mock_update_run_status)
    mock_session = MockSession()

    # test 1: run has cromwell_run_id
    mock_data = MockRunModel(status="queued", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.cancel()
    assert run.status() == "cancelled"

    # test 2: run doesn't have cromwell_run_id
    mock_data = MockRunModel(status="queued", cromwell_run_id=None)
    run = Run(mock_session, mock_data)
    run.cancel()
    assert run.status() == "cancelled"


def test_metadata(monkeypatch):
    class MockMetadata:
        def __init__(self):
            self.data = {"MOCK_METADATA": True}

    def mock_cromwell_get_metadata(self, cromwell_run_id):
        assert type(cromwell_run_id) is str
        mock_metadata = MockMetadata()
        mock_metadata.data["runId"] = cromwell_run_id
        return mock_metadata

    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_metadata", mock_cromwell_get_metadata
    )

    # when no metadata, expect empty dict
    mock_session = MockSession()
    mock_data = MockRunModel(cromwell_run_id=None)
    run = Run(mock_session, mock_data)
    metadata = run.metadata()
    assert metadata is None

    test_cromwell_run_id = "ZZZZ"
    mock_data = MockRunModel(cromwell_run_id=test_cromwell_run_id)
    run = Run(mock_session, mock_data)
    metadata = run.metadata()
    assert metadata.data["runId"] == test_cromwell_run_id


def test_outputs(monkeypatch):
    class MockMetadata:
        def __init__(self):
            self.data = {}

        def outputs(self, relpath=True):
            mock_outputs = {"MOCK-TASK": "MOCK-TASK-OUTPUT"}
            return mock_outputs

    def mock_cromwell_get_metadata(self, cromwell_run_id):
        assert type(cromwell_run_id) is str
        mock_metadata = MockMetadata()
        return mock_metadata

    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_metadata", mock_cromwell_get_metadata
    )

    mock_session = MockSession()

    # an empty dict is expected if Cromwell has not been executed yet
    mock_data = MockRunModel(cromwell_run_id=None)
    run = Run(mock_session, mock_data)
    outputs = run.outputs()
    assert type(outputs) == dict
    assert len(outputs.keys()) == 0

    # if cromwell has run, a dict of outputs is expected
    mock_data = MockRunModel(cromwell_run_id="TEST-CROMWELL-RUN-ID")
    run = Run(mock_session, mock_data)
    outputs = run.outputs()
    assert type(outputs) == dict
    assert outputs["MOCK-TASK"] == "MOCK-TASK-OUTPUT"


def test_outfiles(monkeypatch):
    class MockMetadata:
        def __init__(self):
            self.data = {}

        def outfiles(self, complete=False, relpath=True):
            mock_outfiles = ["./TEST_OUTFILE"]
            return mock_outfiles

    def mock_cromwell_get_metadata(self, cromwell_run_id):
        assert type(cromwell_run_id) is str
        mock_metadata = MockMetadata()
        return mock_metadata

    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_metadata", mock_cromwell_get_metadata
    )

    mock_session = MockSession()

    # an empty list is expected if Cromwell has not been executed yet
    mock_data = MockRunModel(cromwell_run_id=None)
    run = Run(mock_session, mock_data)
    outfiles = run.outfiles()
    assert type(outfiles) == list
    assert len(outfiles) == 0

    # if cromwell has run, a list of files is expected
    mock_data = MockRunModel(cromwell_run_id="TEST-CROMWELL-RUN-ID")
    run = Run(mock_session, mock_data)
    outfiles = run.outfiles()
    assert type(outfiles) == list
    assert outfiles[0] == "./TEST_OUTFILE"


def test_output_manifest(monkeypatch):
    class MockMetadata:
        def __init__(self):
            self.data = {}

        def workflow_root(self):
            return "MOCK-WORKFLOW-ROOT"

        def outfiles(self, complete=False, relpath=True):
            return ["TEST_OUTFILE1", "TEST_OUTFILE2"]

    def mock_cromwell_get_metadata(self, cromwell_run_id):
        assert type(cromwell_run_id) is str
        mock_metadata = MockMetadata()
        return mock_metadata

    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_metadata", mock_cromwell_get_metadata
    )

    mock_session = MockSession()

    # an empty dict is expected if Cromwell has not been executed yet
    mock_data = MockRunModel(cromwell_run_id=None)
    run = Run(mock_session, mock_data)
    output_manifest = run.output_manifest(complete=False)
    assert type(output_manifest) == dict
    assert len(output_manifest.keys()) == 0

    # if cromwell has run, a dict with workflow_root and manifest is expected
    mock_data = MockRunModel(cromwell_run_id="TEST-CROMWELL-RUN-ID")
    run = Run(mock_session, mock_data)
    output_manifest = run.output_manifest(complete=False)
    assert type(output_manifest) == dict
    assert "workflow_root" in output_manifest
    assert "manifest" in output_manifest
    assert type(output_manifest["manifest"]) is list
    assert output_manifest["manifest"][0] == "TEST_OUTFILE1"
    assert output_manifest["workflow_root"] == "MOCK-WORKFLOW-ROOT"


def test_s3_parse_path():
    mock_session = MockSession()
    mock_data = MockRunModel()
    run = Run(mock_session, mock_data)

    test_path = "s3://jaws-bucket/a/b/c/d.txt"
    expected_bucket = "jaws-bucket"
    expected_path = "a/b/c/d.txt"

    bucket, path = run.s3_parse_path(test_path)
    assert bucket == expected_bucket
    assert path == expected_path


def test_inputs(monkeypatch):
    def mock_read_inputs(self):
        example_inputs = {"fasta_file": "CORI/mydata/genome.fasta", "min_score": 95}
        return example_inputs

    monkeypatch.setattr(jaws_site.runs.Run, "read_inputs", mock_read_inputs)

    mock_session = MockSession()
    mock_data = MockRunModel(input_site_id="CORI")
    run = Run(mock_session, mock_data)
    run.config["inputs_dir"] = "/inputs"
    inputs = run.inputs()

    # print(inputs)

    assert inputs["fasta_file"] == "/inputs/CORI/mydata/genome.fasta"
    assert inputs["min_score"] == 95


def test_inputs_fh(monkeypatch):
    def mock_inputs(self):
        return {"KEY": "VALUE"}

    monkeypatch.setattr(jaws_site.runs.Run, "inputs", mock_inputs)

    mock_session = MockSession()
    mock_data = MockRunModel()
    run = Run(mock_session, mock_data)
    fh = run.inputs_fh()
    inputs = json.load(fh)
    assert inputs["KEY"] == "VALUE"


# def test_get_run_inputs(
#    monkeypatch, inputs_files, inputs_files_missing_json, inputs_files_without_zip
# ):
#    def mock__inputs_fh():
#        return None
#
#    def mock__read_file(path):
#        return None
#
#    monkeypatch.setattr(jaws_site.runs.Run, "_inputs_fh", mock__inputs_fh)
#    monkeypatch.setattr(jaws_site.runs.Run, "_read_file", mock__read_file)
#
#    mock_session = MockSession()
#
#    # test 1: valid input
#    data1 = MockRunModel(status="upload complete", submission_id="XXXX")
#    run1 = Run(mock_session, data1)
#    infiles = run1.get_run_inputs()
#    assert infiles[0].read() == "output for XXXX.wdl"
#    assert infiles[1].read() == "output for XXXX.json"
#    assert infiles[2].read().decode() == "output for XXXX.zip"
#    assert infiles[3] is None
#
#    # test 2: invalid input
#    data2 = MockRunModel(status="upload complete", submission_id="YYYY")
#    run2 = Run(mock_session, data2)
#    with pytest.raises(jaws_site.runs.DataError):
#        infiles = run2.get_run_inputs()
#
#    # test 3: valid input, no subworkflows zip
#    data3 = MockRunModel(status="upload complete", submission_id="WWWW")
#    run3 = Run(mock_session, data3)
#    infiles = run3.get_run_inputs()
#    assert infiles[0].read() == "output for WWWW.wdl"
#    assert infiles[1].read() == "output for WWWW.json"
#    assert infiles[2] is None
#    assert infiles[3] is None


def test_submit_run(monkeypatch, inputs_files):
    def mock_cromwell_submit(self, fhs, options):
        return "ABCD-EFGH"

    def mock_get_run_inputs(self):
        return {}, {}

    mock_session = MockSession()
    mock_data = MockRunModel(status="upload complete")
    run = Run(mock_session, mock_data)

    monkeypatch.setattr(jaws_site.cromwell.Cromwell, "submit", mock_cromwell_submit)
    monkeypatch.setattr(Run, "get_run_inputs", mock_get_run_inputs)
    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)

    run.submit_run()


@pytest.fixture
def mock_path(tmp_path):
    def tmp_path_mock(*args, **kwargs):
        return (tmp_path / "cwd/inputs/jaws/XXXX").as_posix()

    return tmp_path_mock


def test_check_run_cromwell_status(monkeypatch):

    monkeypatch.setattr(Run, "_update_run_status", mock__update_run_status)
    monkeypatch.setattr(Run, "_insert_run_log", mock__insert_run_log)

    def mock_get_status_running(self, run_id):
        return "Running"

    def mock_get_status_succeeded(self, run_id):
        return "Succeeded"

    def mock_get_status_failed(self, run_id):
        return "Failed"

    def mock_did_run_start_true(self):
        return True

    def mock_did_run_start_false(self):
        return False

    mock_session = MockSession()

    # test: submitted -> queued
    mock_data = MockRunModel(status="queued")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.runs.Run, "did_run_start", mock_did_run_start_false)
    run.check_run_cromwell_status()
    assert run.data.status == "queued"

    # test: queued -> queued
    mock_data = MockRunModel(status="queued")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    run.check_run_cromwell_status()
    assert run.data.status == "queued"

    # test: queued -> running
    mock_data = MockRunModel(status="queued")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.runs.Run, "did_run_start", mock_did_run_start_true)
    run.check_run_cromwell_status()
    assert run.data.status == "running"

    # test: queued -> succeeded
    mock_data = MockRunModel(status="queued")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    run.check_run_cromwell_status()
    assert run.data.status == "succeeded"

    # test: queued -> failed
    mock_data = MockRunModel(status="queued")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    run.check_run_cromwell_status()
    assert run.data.status == "failed"

    # test: running -> succeeded
    mock_data = MockRunModel(status="running")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    run.check_run_cromwell_status()
    assert run.data.status == "succeeded"

    # test: running -> failed
    mock_data = MockRunModel(status="running")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    run.check_run_cromwell_status()
    assert run.data.status == "failed"


def test_run_log(monkeypatch):
    def mock_select_rows(self):
        mock_results = [
            [123, "queued", "running", "2022-05-03 11:06:05", None, True],
            [123, "running", "failed", "2022-05-03 11:18:44", "failure reason", True],
        ]
        return mock_results

    monkeypatch.setattr(RunLog, "_select_rows", mock_select_rows)

    exp_results = [
        {
            "status_from": "queued",
            "status_to": "running",
            "timestamp": "2022-05-03 11:06:05",
            "reason": None,
        },
        {
            "status_from": "running",
            "status_to": "failed",
            "timestamp": "2022-05-03 11:18:44",
            "reason": "failure reason",
        },
    ]

    mock_session = None
    run_id = 123
    run_log = RunLog(mock_session, run_id)
    obs_results = run_log.logs()

    assert bool(DeepDiff(obs_results, exp_results, ignore_order=True)) is False


@pytest.mark.skip(reason="getting rid of singleton config.conf")
def test_summary():
    mock_session = MockSession()
    mock_data = MockRunModel(
        id=123,
        user_id="jdoe",
        email="johndoe@lbl.gov",
        submitted=this_date,
        updated=this_date,
        status="download complete",
        result="succeeded",
    )
    run = Run(mock_session, mock_data)
    exp_results = {
        "run_id": 123,
        "user_id": "jdoe",
        "email": "johndoe@lbl.gov",
        "submitted": this_date.strftime("%Y-%m-%d %H:%M:%S"),
        "updated": this_date.strftime("%Y-%m-%d %H:%M:%S"),
        "status": "download complete",
        "result": "succeeded",
        "compute_site_id": "EAGLE",
        "status_detail": "",
    }
    run = Run(mock_session, mock_data)
    obs_results = run.summary()
    assert bool(DeepDiff(obs_results, exp_results, ignore_order=True)) is False


def test_from_params(mock_sqlalchemy_session):
    params = {
        "user_id": "testuser",
        "submission_id": 123,
        "max_ram_gb": 11,
        "caching": True,
        "input_site_id": "JGI",
        "compute_site_id": "JGI",
        "wdl_file": "test.wdl",
        "json_file": "test.json",
        "tag": "atag",
        "webhook": "xxx",
        "manifest": ["test", "test2"],
        "manifest_json": "test.json",
        "run_id": 123,
    }
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    run.from_params(mock_sqlalchemy_session, params)

    assert mock_sqlalchemy_session.data.session["add"] is True
    assert mock_sqlalchemy_session.data.session["commit"] is True
    mock_sqlalchemy_session.clear()

    # Test add exception
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    with pytest.raises(RunDbError):
        run.from_params(mock_sqlalchemy_session, params)

    assert mock_sqlalchemy_session.data.session["rollback"] is True
    mock_sqlalchemy_session.clear()

    # Test general Exception
    mock_sqlalchemy_session.data.raise_exception = True
    with pytest.raises(Exception):
        run.from_params(mock_sqlalchemy_session, params)

    assert mock_sqlalchemy_session.data.session["rollback"] is True
    mock_sqlalchemy_session.clear()

    # Test invalid data types for run model
    params = {
        "run_id": 123,
        "user_id": 123,
        "submission_id": "123",
        "caching": 123,
    }
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    with pytest.raises(RunDbError):
        run.from_params(mock_sqlalchemy_session, params)


def test_from_id(mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    run.from_id(mock_sqlalchemy_session, 123)
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test raising Exception
    mock_sqlalchemy_session.data.raise_exception = True
    with pytest.raises(Exception):
        run.from_id(mock_sqlalchemy_session, 123)
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test raising SQLAlchemyError
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    with pytest.raises(RunDbError):
        run.from_id(mock_sqlalchemy_session, 123)
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()


def test_from_cromwell_run_id(mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    run.from_cromwell_run_id(mock_sqlalchemy_session, "myid")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test raising Exception
    mock_sqlalchemy_session.data.raise_exception = True
    with pytest.raises(Exception):
        run.from_cromwell_run_id(mock_sqlalchemy_session, "myid")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test raising SQLAlchemyError
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    with pytest.raises(RunDbError):
        run.from_cromwell_run_id(mock_sqlalchemy_session, "myid")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test RunNotFound
    mock_sqlalchemy_session.data.raise_exception_noresultfound = True
    with pytest.raises(RunNotFoundError):
        run.from_cromwell_run_id(mock_sqlalchemy_session, "myid")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()


def test_summary2(mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.summary()
    assert ret["run_id"] == "99"
    assert ret["user_id"] == "test_user"
    assert ret["cromwell_run_id"] == "myid"
    assert ret["result"] == "succeeded"


def test_did_run_start(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.did_run_start()
    assert ret is True


def test_task_summary(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.task_summary()
    assert ret == [
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
        }
    ]

    data = initRunModel(status="succeeded")
    data.cromwell_run_id = None
    run = Run(mock_sqlalchemy_session, data)
    ret = run.task_summary()
    assert len(ret) == 0


def test_task_log(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.task_log()
    assert ret == {
        "name": "test",
        "cached": False,
        "execution_status": "succeeded",
        "queue_start": "01-01-2022",
        "run_start": "01-01-2022",
        "run_end": "01-01-2022",
        "queue_duration": "01-01-2022",
        "run_duration": "01-01-2022",
    }

    data = initRunModel(status="succeeded")
    data.cromwell_run_id = None
    run = Run(mock_sqlalchemy_session, data)
    ret = run.task_log()
    assert len(ret) == 0


def test_report(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.report()
    ret == {
        "run_id": "99",
        "user_id": "test_user",
        "cromwell_run_id": "myid",
        "submitted": "2022-09-15 00:06:13",
        "updated": "2022-09-15 00:06:13",
        "status": "succeeded",
        "result": "succeeded",
        "workflow_name": "unknown",
        "site_id": "EAGLE",
        "tasks": [
            {
                "name": "test",
                "shard_index": "shard_index",
                "attempt": "attempt",
                "cached": "cached",
                "execution_status": "succeeded",
                "result": "result",
                "failure_message": "failure_message",
                "queue_start": "01-01-2022",
                "run_start": "01-01-2022",
                "run_end": "01-01-2022",
                "call_root": "call_root",
                "requested_time": "01-01-2022",
                "requested_cpu": 15,
                "requested_memory": "100",
                "cromwell_job_id": "123",
                "queuetime_sec": 0,
                "runtime_sec": 0,
                "walltime_sec": 0,
            }
        ],
    }


def test_cancel2(mock_metadata, mock_sqlalchemy_session, monkeypatch):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.cancel()
    assert ret is None

    def mock_cromwell_abort(self, cid):
        pass

    monkeypatch.setattr(Cromwell, "abort", mock_cromwell_abort)

    data = initRunModel(status="submitted")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.cancel()
    assert ret is None

    # Test exception in cromwell.abort
    def mock_cromwell_abort(self, cid):
        raise CromwellError

    monkeypatch.setattr(Cromwell, "abort", mock_cromwell_abort)

    data = initRunModel(status="submitted")
    run = Run(mock_sqlalchemy_session, data)
    with pytest.raises(CromwellError):
        run.cancel()

    # Test RunDbError
    def mock_update_run_status(self, status):
        raise Exception

    monkeypatch.setattr(Run, "update_run_status", mock_update_run_status)

    with pytest.raises(RunDbError):
        run.cancel()


def test_workflow_root(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.workflow_root()
    assert ret == "/root"

    data.cromwell_run_id = None
    run = Run(mock_sqlalchemy_session, data)
    ret = run.workflow_root()
    assert ret is None


def test_errors(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.errors()
    assert ret == "error"

    data.cromwell_run_id = None
    run = Run(mock_sqlalchemy_session, data)
    ret = run.errors()
    assert ret == {}


def test_running_tasks(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.running_tasks()
    assert ret == "running"

    data.cromwell_run_id = None
    run = Run(mock_sqlalchemy_session, data)
    ret = run.running_tasks()
    assert ret == {}


def test__read_file_nfs(
    config_file,
    file_not_found_config,
    config_file_zero,
    mock_metadata,
    mock_sqlalchemy_session,
):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run._read_file_nfs(config_file)
    assert isinstance(ret, io.StringIO)

    # Test OSError
    with pytest.raises(OSError):
        run._read_file_nfs(file_not_found_config)

    # Test zero file
    with pytest.raises(OSError):
        run._read_file_nfs(config_file_zero)


def test__read_file(config_file, monkeypatch, mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)

    def mock__read_file_s3(self, path, binary):
        return "s3"

    monkeypatch.setattr(Run, "_read_file_s3", mock__read_file_s3)

    def mock__read_file_nfs(self, path, binary):
        return "nfs"

    monkeypatch.setattr(Run, "_read_file_nfs", mock__read_file_nfs)

    ret = run._read_file(config_file)
    assert ret == "nfs"

    run.config["inputs_dir"] = "s3://"
    ret = run._read_file(config_file)
    assert ret == "s3"


def test_read_inputs(monkeypatch, mock_sqlalchemy_session):
    import os

    tests_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = "inputs.json"

    def mock__read_file(self, path):
        fname = f"{tests_dir}/{input_file}"
        fh = open(fname, "r")
        return fh

    monkeypatch.setattr(Run, "_read_file", mock__read_file)

    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.read_inputs()
    assert ret == {
        "main_wdl.fastq": "/global/JAWS/data/veryshort.fastq.bz2",
        "main_wdl.reference": "/global/JAWS/data/DOE_UTEX.polished.t635masked.fasta",
        "main_wdl.bbmap_shard_wf.chunk_size": "100000000",
        "main_wdl.bbtools_mem": "10G",
    }


def test_get_run_inputs(monkeypatch, mock_sqlalchemy_session):
    import os

    tests_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = "inputs.json"

    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)

    # Test RunFileNotFoundError
    with pytest.raises(RunFileNotFoundError):
        run.get_run_inputs()

    def mock_inputs_fh(self):
        fname = f"{tests_dir}/{input_file}"
        fh = open(fname, "r")
        return fh

    monkeypatch.setattr(Run, "inputs_fh", mock_inputs_fh)

    with pytest.raises(RunFileNotFoundError):
        run.get_run_inputs()
