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
    RunInputError,
    send_run_status_logs,
)
from datetime import datetime
from tests.conftest import MockSession, MockRunModel, initRunModel
from jaws_site.cromwell import Cromwell, CromwellError, CromwellServiceError
from jaws_rpc.rpc_client_basic import RpcClientBasic
from unittest.mock import patch

# from unittest.mock import MagicMock
import io


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


@pytest.fixture
def mock_rpc_client():
    sent_requests = []

    def record_sent_request(*args):
        sent_requests.append({"args": args})
        return {"success": True}

    with patch("jaws_rpc.rpc_client_basic.RpcClientBasic") as mock_client:
        mock_client.return_value.request.side_effect = record_sent_request
        yield mock_client.return_value, sent_requests


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
    assert run.data.status == "cancel"

    # test2: run was already marked to be cancelled
    mock_data = MockRunModel(status="cancel", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.data.status == "cancel"

    # test3: run was already cancelled
    mock_data = MockRunModel(status="cancelled", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.data.status == "cancelled"

    # test 4: run is active and has cromwell_run_id
    mock_data = MockRunModel(status="queued", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.data.status == "cancel"

    # test 5: run is succeeded
    mock_data = MockRunModel(status="succeeded", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    with pytest.raises(jaws_site.runs.RunInputError):
        run.mark_to_cancel()

    # test 6: run is ready to submit to Cromwell
    mock_data = MockRunModel(status="ready", cromwell_run_id=None)
    run = Run(mock_session, mock_data)
    run.mark_to_cancel()
    assert run.data.status == "cancelled"


def test_cancel(monkeypatch):
    def mock_cromwell_abort(self, cromwell_run_id):
        assert type(cromwell_run_id) is str
        return {"id": cromwell_run_id, "status": "Aborting"}

    def mock_update_run_status(self, new_status):
        assert new_status is not None
        assert type(new_status) == str
        assert new_status != self.data.status
        self.data.status = new_status

    monkeypatch.setattr(jaws_site.cromwell.Cromwell, "abort", mock_cromwell_abort)
    monkeypatch.setattr(Run, "update_run_status", mock_update_run_status)
    mock_session = MockSession()

    mock_data = MockRunModel(status="cancel", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.cancel()
    assert run.data.status == "cancelled"

    # test 2: cromwell is unavailable
    def mock_cromwell_abort_raises(self, cromwell_run_id):
        raise CromwellServiceError("Cromwell service in unavailable")

    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "abort", mock_cromwell_abort_raises
    )

    mock_data = MockRunModel(status="cancel", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.cancel()
    assert run.data.status == "cancel"


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


def test_submit_run(mock_sqlalchemy_session, monkeypatch, inputs_files):
    def mock_cromwell_submit(self, fhs, options):
        return "ABCD-EFGH"

    def mock_get_run_inputs(self):
        return {}, {}

    def mock_max_active_runs_exceeded(*args):
        return False

    mock_data = MockRunModel(status="upload complete")
    run = Run(mock_sqlalchemy_session, mock_data)

    monkeypatch.setattr(jaws_site.cromwell.Cromwell, "submit", mock_cromwell_submit)
    monkeypatch.setattr(Run, "get_run_inputs", mock_get_run_inputs)
    monkeypatch.setattr(
        jaws_site.runs, "max_active_runs_exceeded", mock_max_active_runs_exceeded
    )

    run.submit_run()
    assert mock_sqlalchemy_session.data.session["commit"] is True
    assert mock_sqlalchemy_session.data.session["close"] is False


def test_resubmit_run(mock_sqlalchemy_session):
    # test 1: fail active run
    mock_data = MockRunModel(status="running")
    run = Run(mock_sqlalchemy_session, mock_data)
    with pytest.raises(RunInputError):
        run.resubmit()

    # test 2: success for finished run
    mock_data = MockRunModel(status="finished")
    run = Run(mock_sqlalchemy_session, mock_data)
    run.resubmit()
    assert mock_sqlalchemy_session.data.session["commit"] is True
    assert mock_sqlalchemy_session.data.session["close"] is False


def test_max_active_runs_exceeded(mock_sqlalchemy_session):
    # test user has active run but under exceeded threshold
    mock_sqlalchemy_session.output([{"run_id": 123}, {"run_id": 456}])
    obs_result = jaws_site.runs.max_active_runs_exceeded(
        mock_sqlalchemy_session, 123, 10
    )
    assert obs_result is False

    # test user has active run and exceeds threshold
    obs_result = jaws_site.runs.max_active_runs_exceeded(
        mock_sqlalchemy_session, 123, 1
    )
    assert obs_result is True

    # test user does not have any active run
    mock_sqlalchemy_session.clear()
    obs_result = jaws_site.runs.max_active_runs_exceeded(
        mock_sqlalchemy_session, 123, 1
    )
    assert obs_result is False


@pytest.fixture
def mock_path(tmp_path):
    def tmp_path_mock(*args, **kwargs):
        return (tmp_path / "cwd/inputs/jaws/XXXX").as_posix()

    return tmp_path_mock


def test_check_cromwell_run_status(monkeypatch, mock_metadata):
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
    mock_data = MockRunModel(status="queued", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.runs.Run, "did_run_start", mock_did_run_start_false)
    run.check_cromwell_run_status()
    assert run.data.status == "queued"

    # test: queued -> queued
    mock_data = MockRunModel(status="queued", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    run.check_cromwell_run_status()
    assert run.data.status == "queued"

    # test: queued -> running
    mock_data = MockRunModel(status="queued", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_running
    )
    monkeypatch.setattr(jaws_site.runs.Run, "did_run_start", mock_did_run_start_true)
    run.check_cromwell_run_status()
    assert run.data.status == "running"

    # test: queued -> succeeded
    mock_data = MockRunModel(status="queued", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    run.check_cromwell_run_status()
    assert run.data.status == "succeeded"

    # test: queued -> failed
    mock_data = MockRunModel(status="queued", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    run.check_cromwell_run_status()
    assert run.data.status == "failed"

    # test: running -> succeeded
    mock_data = MockRunModel(status="running", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_succeeded
    )
    run.check_cromwell_run_status()
    assert run.data.status == "succeeded"

    # test: running -> failed
    mock_data = MockRunModel(status="running", workflow_root="/x")
    run = Run(mock_session, mock_data)
    monkeypatch.setattr(
        jaws_site.cromwell.Cromwell, "get_status", mock_get_status_failed
    )
    run.check_cromwell_run_status()
    assert run.data.status == "failed"

    # test get metadata, workflow_root
    mock_data = MockRunModel(status="submitted", cromwell_run_id="ABCD")
    run = Run(mock_session, mock_data)
    run.check_cromwell_run_status()
    assert run.data.status == "queued"
    assert run.data.workflow_name == "unknown"
    assert run.data.workflow_root == "/data/cromwell-executions/example/ABCD"


def test_run_log(monkeypatch):
    def mock_select_rows(self):
        mock_results = [
            [3825, 123, "queued", "running", "2022-05-03 11:06:05", None, True],
            [
                3895,
                123,
                "running",
                "failed",
                "2022-05-03 11:18:44",
                "failure reason",
                True,
            ],
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
    assert mock_sqlalchemy_session.data.session["close"] is False
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
    run.from_cromwell_run_id(mock_sqlalchemy_session, "ABCD-EFGH")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test raising Exception
    mock_sqlalchemy_session.data.raise_exception = True
    with pytest.raises(Exception):
        run.from_cromwell_run_id(mock_sqlalchemy_session, "ABCD-EFGH")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test raising SQLAlchemyError
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    with pytest.raises(RunDbError):
        run.from_cromwell_run_id(mock_sqlalchemy_session, "ABCD-EFGH")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test RunNotFound
    mock_sqlalchemy_session.data.raise_exception_noresultfound = True
    with pytest.raises(RunNotFoundError):
        run.from_cromwell_run_id(mock_sqlalchemy_session, "ABCD-EFGH")
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()


def test_summary2(mock_sqlalchemy_session):
    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.summary()
    assert ret["run_id"] == "99"
    assert ret["user_id"] == "test_user"
    assert ret["cromwell_run_id"] == "ABCD-EFGH"
    assert ret["result"] == "succeeded"


def test_did_run_start(mock_metadata, mock_sqlalchemy_session):

    data = initRunModel(cromwell_run_id=None)
    run = Run(mock_sqlalchemy_session, data)
    ret = run.did_run_start()
    assert ret is False


def test_task_log(mock_metadata, mock_sqlalchemy_session):
    data = initRunModel(cromwell_run_id=None)
    run = Run(mock_sqlalchemy_session, data)
    ret = run.task_log()
    assert len(ret) == 0


def test_summary(mock_metadata, mock_sqlalchemy_session, monkeypatch):

    data = initRunModel(status="succeeded")
    run = Run(mock_sqlalchemy_session, data)
    ret = run.summary()
    ret == {
        "run_id": "99",
        "user_id": "test_user",
        "cromwell_run_id": "ABCD-EFGH",
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


def test_publish_report(mock_metadata, mock_sqlalchemy_session, monkeypatch):
    def mock__read_json_file(self, path):
        contents = {
            "run_id": "99",
        }
        return contents

    def mock_rpc_client_request(self, payload):
        successful_response = {"result": {}}
        return successful_response

    def mock_update_run_status(self, new_status):
        assert new_status

    monkeypatch.setattr(Run, "_read_json_file", mock__read_json_file)
    monkeypatch.setattr(RpcClientBasic, "request", mock_rpc_client_request)

    data = initRunModel(status="complete")
    run = Run(mock_sqlalchemy_session, data)
    run.publish_report()


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


def test_send_run_status_logs(mock_sqlalchemy_session, mock_rpc_client, tmpdir):
    mock_rpc_client, sent_requests = mock_rpc_client
    cromwell_run_id = "6649faab-e485-48ed-8b23-bc25da8b509e"
    workflow_root = tmpdir.mkdir(cromwell_run_id)
    outputs_manifest_contents = [
        "metadata.json",
        "errors.json",
        "outputs.json",
        "outmanifest.json",
        "task_summary.json",
    ]
    outputs_manifest_file = workflow_root.join("output_manifest.json")
    with open(outputs_manifest_file, "w") as f:
        json.dump(outputs_manifest_contents, f)

    ready_time = datetime.strptime("2023-06-10 12:00:00", "%Y-%m-%d %H:%M:%S")
    running_time = datetime.strptime("2023-06-10 12:10:00", "%Y-%m-%d %H:%M:%S")
    completion_time = datetime.strptime("2023-06-10 12:20:00", "%Y-%m-%d %H:%M:%S")
    log_query_results = [
        {
            "sent": False,
            "run_id": 1,
            "status_from": "ready",
            "status_to": "submitted",
            "timestamp": ready_time,
            "reason": "Test reason 1",
        },
        {
            "sent": False,
            "run_id": 1,
            "status_from": "submitted",
            "status_to": "queued",
            "timestamp": running_time,
            "reason": "Test reason 2",
        },
        {
            "sent": False,
            "run_id": 1,
            "status_from": "running",
            "status_to": "complete",
            "timestamp": completion_time,
            "reason": "Test reason 3",
        },
    ]
    run_query_results = [
        {
            "id": 1,
            "run_id": 1,
            "cromwell_run_id": cromwell_run_id,
            "workflow_name": "fq_count",
            "workflow_root": str(workflow_root),
        }
    ]
    mock_sqlalchemy_session.output(log_query_results)
    mock_sqlalchemy_session.output(run_query_results)

    send_run_status_logs(mock_sqlalchemy_session, mock_rpc_client)

    assert len(sent_requests) == 3
    request_bodies = [sent_request["args"][1] for sent_request in sent_requests]
    assert request_bodies[0]["cromwell_run_id"] == cromwell_run_id
    assert request_bodies[1]["workflow_name"] == "fq_count"
    assert request_bodies[1]["workflow_root"] == str(workflow_root)
    assert request_bodies[2]["output_manifest"] == [
        "metadata.json",
        "errors.json",
        "outputs.json",
        "outmanifest.json",
        "task_summary.json",
    ]


def test_task_summary(
    requests_mock, mock_metadata, mock_sqlalchemy_session, monkeypatch
):
    def mock_task_log(self):
        return {
            "header": [
                "TASK_DIR",
                "STATUS",
                "QUEUE_START",
                "RUN_START",
                "RUN_END",
                "RC",
                "QUEUE_DUR",
                "RUN_DUR",
            ],
            "data": [
                [
                    "call-test",
                    "done",
                    "01-01-2022 01:00:00",
                    "01-01-2022 01:01:00",
                    "01-01-2022 01:11:00",
                    0,
                    "00:01:00",
                    "00:10:00",
                ]
            ],
        }

    monkeypatch.setattr(Run, "task_log", mock_task_log)

    data = initRunModel()
    run = Run(mock_sqlalchemy_session, data)

    actual = run.task_summary()
    expected = [
        {
            "name": "test",
            "shard_index": "-1",
            "attempt": 1,
            "cached": False,
            "job_id": "123",
            "execution_status": "done",
            "result": "succeeded",
            "failure_message": None,
            "status": "done",
            "queue_start": "01-01-2022 01:00:00",
            "run_start": "01-01-2022 01:01:00",
            "run_end": "01-01-2022 01:11:00",
            "rc": 0,
            "queue_duration": "00:01:00",
            "run_duration": "00:10:00",
            "call_root": "/scratch/cromwell-executions/testWorkflow/ABCD-EFGH/call-test",
            "requested_time": "00:30:00",
            "requested_cpu": 15,
            "requested_memory": "10 GB",
        },
    ]
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False
