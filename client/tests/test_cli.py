import pytest
import os
import requests
import shutil
import subprocess
import json
import click.testing
from jaws_client import cli


# flake8: noqa
HISTORY = [
        {
        "id": "33",
        "input_site_id": "CORI",
        "json_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.json",
        "result": "succeeded",
        "site_id": "CORI",
        "status": "download complete",
        "status_detail": "The run output (succeeded or failed) has been returned to the user.",
        "submitted": "2021-01-01 11:00:00",
        "tag": "none",
        "updated": "2021-01-01 12:00:00",
        "wdl_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.wdl"
        },
        {
        "id": "34",
        "input_site_id": "CORI",
        "json_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.json",
        "result": "succeeded",
        "site_id": "CORI",
        "status": "download complete",
        "status_detail": "The run output (succeeded or failed) has been returned to the user.",
        "submitted": "2021-07-13 14:00:00",
        "tag": "none",
        "updated": "2021-07-13 14:51:55",
        "wdl_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.wdl"
        }
]


QUEUE = [
        {
        "id": "33",
        "input_site_id": "CORI",
        "json_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.json",
        "result": "succeeded",
        "site_id": "CORI",
        "status": "download complete",
        "status_detail": "The run output (succeeded or failed) has been returned to the user.",
        "submitted": "2021-01-01 11:00:00",
        "tag": "none",
        "updated": "2021-01-01 12:00:00",
        "wdl_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.wdl"
        },
        {
        "id": "34",
        "input_site_id": "CORI",
        "json_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.json",
        "result": "succeeded",
        "site_id": "CORI",
        "status": "download complete",
        "status_detail": "The run output (succeeded or failed) has been returned to the user.",
        "submitted": "2021-07-13 14:00:00",
        "tag": "none",
        "updated": "2021-07-13 14:51:55",
        "wdl_file": "/global/cscratch1/sd/jaws/jfroula/jaws-health-test/fq_count.wdl"
        }
]


WORKFLOW_METADATA = {
    "run_id": "36",
    "actualWorkflowLanguage": "WDL",
    "actualWorkflowLanguageVersion": "draft-2",
    "calls": {
        "bbtools.alignment": [
            {
                "attempt": 1,
                "backend": "JTM",
                "backendStatus": "Running",
                "callCaching": {
                    "allowResultReuse": "true",
                    "effectiveCallCachingMode": "ReadAndWriteCache",
                    "hit": "false",
                    "result": "Cache Miss",
                },
                "callRoot": "/tmp/jaws-dev/cromwell/cromwell-executions/bbtools/f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e/call-alignment",  # noqa
                "commandLine": "shifterimg pull jfroula/bbtools:1.2.1 && \\\nshifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=/tmp/jaws-dev/cromwell/cromwell-executions/bbtools/f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e/call-alignment/inputs/-1159340168/5min_reads.fq ref=/tmp/jaws-dev/cromwell/cromwell-executions/bbtools/f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e/call-alignment/inputs/-1159340168/5min_ref.fasta out=test.sam",  # noqa
                "executionStatus": "Running",
                "inputs": {
                    "fasta": "/global/dna/shared/data/jfroula/JAWS/data/5min_ref.fasta",
                    "fastq": "/global/dna/shared/data/jfroula/JAWS/data/5min_reads.fq",
                },
                "jobId": "778690",
                "runtimeAttributes": {
                    "account": "fungalp",
                    "cluster": "cori",
                    "constraint": "haswell",
                    "continueOnReturnCode": "0",
                    "cpu": "1",
                    "failOnStderr": "false",
                    "maxRetries": "0",
                    "mem": "0G",
                    "node": "1",
                    "nwpn": "1",
                    "poolname": "small",
                    "qos": "genepool",
                    "shared": "1",
                    "time": "00:00:00",
                },
                "shardIndex": -1,
                "start": "2020-04-03T20:32:52.938Z",
                "stderr": "/tmp/jaws-dev/cromwell/cromwell-executions/bbtools/f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e/call-alignment/execution/stderr",  # noqa
                "stdout": "/tmp/jaws-dev/cromwell/cromwell-executions/bbtools/f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e/call-alignment/execution/stdout",  # noqa
            }
        ]
    },
    "id": "f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e",
    "inputs": {
        "bbtools.reads": "/global/dna/shared/data/jfroula/JAWS/data/5min_reads.fq",
        "bbtools.ref": "/global/dna/shared/data/jfroula/JAWS/data/5min_ref.fasta",
    },
    "labels": {
        "cromwell-workflow-id": "cromwell-f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e",
        "username": "mamelara",
    },
    "outputs": {},
    "start": "2020-04-03T20:32:50.080Z",
    "status": "Running",
    "submission": "2020-04-03T20:32:49.265Z",
    "submittedFiles": {
        "inputs": '{"bbtools.reads":"/global/dna/shared/data/jfroula/JAWS/data/5min_reads.fq","bbtools.ref":"/global/dna/shared/data/jfroula/JAWS/data/5min_ref.fasta"}',
        "labels": '{"username": "mamelara"}',
        "options": "{\n\n}",
        "root": "",
        "workflow": 'workflow bbtools { File reads\n    File ref\n\n    call alignment {\n       input: fastq=reads,\n              fasta=ref\n    }\n    call samtools {\n       input: sam=alignment.sam\n   }\n}\n\ntask alignment {\n    File fastq\n    File fasta\n\n    command <<<\n        shifterimg pull jfroula/bbtools:1.2.1 && \\\n        shifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=${fastq} ref=${fasta} out=test.sam\n    >>>\n    output {\n       File sam = "test.sam"\n    }\n}\n\ntask samtools {\n    File sam\n\n    command {\n       shifter --image=jfroula/bbtools:1.2.1 samtools view -b -F0x4 ${sam} | shifter --image=jfroula/bbtools:1.2.1 samtools sort - > test.sorted.bam\n    }\n    output {\n       File bam = "test.sorted.bam"\n    }\n\t#runtime {\n\t#  docker: "jfroula/bbtools:1.2.1"\n\t#}\n}\n',
        "workflowUrl": "",
    },
    "workflowName": "bbtools",
    "workflowProcessingEvents": [
        {
            "cromwellId": "cromid-22ffeb7",
            "cromwellVersion": "45.1",
            "description": "PickedUp",
            "timestamp": "2020-04-03T20:32:50.048Z",
        }
    ],
    "workflowRoot": "/tmp/jaws-dev/cromwell/cromwell-executions/bbtools/f256dc24-971f-4ba2-9e1f-6f7a53ec3e9e",
}


RUN_LOG_JSON = [
    [
        "created",
        "uploading",
        "2020-06-08 06:28:36",
        "upload_task_id=4884e9a8-a951-11ea-9a3b-0255d23c44ef",
    ],
    ["uploading", "upload complete", "2020-06-08 06:28:50", ""],
    [
        "upload complete",
        "submitted",
        "2020-06-08 06:29:01",
        "cromwell_run_id=5d1ba0bd-ef40-42dd-b33b-0c4c31174b76",
    ],
    ["submitted", "running", "2020-06-08 06:29:11", ""],
    ["running", "succeeded", "2020-06-08 06:29:21", ""],
    ["succeeded", "ready", "2020-06-08 06:29:35", ""],
    [
        "ready",
        "downloading",
        "2020-06-08 06:29:47",
        "Globus download_task_id=72a56d84-a951-11ea-bee5-0e716405a293",
    ],
    ["downloading", "download complete", "2020-06-08 06:30:34", ""],
    ["download complete", "finished", "2020-06-08 06:30:44", ""],
]

RUN_LOG_TEXT = (
    "#STATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON\n"
    "created\tuploading\t2020-06-08 06:28:36\tupload_task_id=4884e9a8-a951-11ea-9a3b-0255d23c44ef\n"
    "uploading\tupload complete\t2020-06-08 06:28:50\t\n"
    "upload complete\tsubmitted\t2020-06-08 06:29:01\tcromwell_run_id=5d1ba0bd-ef40-42dd-b33b-0c4c31174b76\n"
    "submitted\trunning\t2020-06-08 06:29:11\t\n"
    "running\tsucceeded\t2020-06-08 06:29:21\t\n"
    "succeeded\tready\t2020-06-08 06:29:35\t\n"
    "ready\tdownloading\t2020-06-08 06:29:47\tGlobus download_task_id=72a56d84-a951-11ea-bee5-0e716405a293\n"
    "downloading\tdownload complete\t2020-06-08 06:30:34\t\n"
    "download complete\tfinished\t2020-06-08 06:30:44\t\n"
)

TASK_LOG_JSON = [
    [
        "runblastplus_sub.task1",
        1,
        43,
        "ready",
        "queued",
        "2020-06-10 13:42:44",
        "",
        "The job was received by JTM-manager and sent to JTM-worker",
    ],
    [
        "runblastplus_sub.task2",
        1,
        44,
        "queued",
        "pending",
        "2020-06-10 13:43:36",
        "",
        "The job was receive by JTM-worker and is awaiting resources",
    ],
]

TASK_LOG_TEXT = (
    "#TASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON\tSTATUS_DETAIL\n"
    "runblastplus_sub.task1\t1\t43\tready\tqueued\t2020-06-10 13:42:44\t\tThe job was received by JTM-manager and sent to JTM-worker\n"
    "runblastplus_sub.task2\t1\t44\tqueued\tpending\t2020-06-10 13:43:36\t\tThe job was receive by JTM-worker and is awaiting resources\n"
)

SUBMISSION = {
    "output_dir": "/global/homes/m/mamelara/out",
    "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    "run_id": 35,
    "site_id": "NERSC",
    "status": "uploading",
    "submission_id": "65f2f4df-2a6c-4881-a3b0-3141107ac668",
    "upload_task_id": "1aba54ac-7695-11ea-9615-0afc9e7dd773",
}


TASK_STATUS_JSON = [
    ["bbtools.alignment", 1, 432, "Queued", "Running", "2020-04-03 20:32:52", ""]
]

TASK_STATUS_TEXT = (
    "#TASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON\n"
    "bbtools.alignment\t1\t432\tqueued\trunning\t2020-04-03 20:32:52\t\n"
)


class MockResponse:
    def __init__(self, result, status_code):
        self.result = result
        self.status_code = status_code

    def json(self):
        return self.result

    @property
    def text(self):
        return str(self.result)


def test_cli_queue(monkeypatch, configuration):
    def post_queue(url, data={}, files={}, headers=None):
        return MockResponse(QUEUE, 200)

    monkeypatch.setattr(requests, "post", post_queue)

    runner = click.testing.CliRunner()
    result = runner.invoke(cli.main, ["queue"])
    assert result.exit_code == 0
    for task_id in ["33", "34"]:
        assert task_id in result.output

    # this checks that there is 8hrs difference when we are not in daylight savings (nov 8 - march 13)
    # utc: 2021-01-01 11:00:00
    assert('2021-01-01 03:00:00' in result.output)

    # this checks that there is 7hrs difference when we are in daylight savings(march 14 - nov 7)
    # utc: 2021-07-13 14:00:00
    assert('2021-07-13 07:00:00' in result.output)


def test_cli_history(monkeypatch, configuration):
    def post_history(url, data={}, files={}, headers=None):
        return MockResponse(HISTORY, 200)

    monkeypatch.setattr(requests, "post", post_history)
    runner = click.testing.CliRunner()
    result = runner.invoke(cli.main, ["history"])
    assert result.exit_code == 0

    for task_id in ["33", "34"]:
        assert task_id in result.output

    # this checks that there is 8hrs difference when we are not in daylight savings (nov 8 - march 13)
    # utc: 2021-01-01 11:00:00
    assert('2021-01-01 03:00:00' in result.output)

    # this checks that there is 7hrs difference when we are in daylight savings(march 14 - nov 7)
    # utc: 2021-07-13 14:00:00
    assert('2021-07-13 07:00:00' in result.output)


def test_cli_status(monkeypatch, configuration):
    def mock_status_get(url, headers={}):
        """
        Returns a list of attributes from a run ordered as the following:
        run id, submission date, submission id, upload id
        """
        job_status = {"status": "Running",'submitted': '2021-07-06 22:40:04', 'tag': None, 'updated': '2021-07-06 22:42:55'}
        if url.endswith('/complete'):
            job_status["output_dir"] = "/foo/bar"
        return MockResponse(job_status, 200)

    monkeypatch.setattr(requests, "get", mock_status_get)
    runner = click.testing.CliRunner()

    result = runner.invoke(cli.main, ["status", "36", "--verbose"])
    assert result.exit_code == 0
    assert "Running" in result.output
    assert "/foo/bar" in result.output

    result = runner.invoke(cli.main, ["status", "36"])
    assert result.exit_code == 0
    assert "Running" in result.output
    assert "/foo/bar" not in result.output

    assert('2021-07-06 15:40:04' in result.output or '2021-07-06 14:40:04' in result.output)


def test_cli_metadata(monkeypatch, configuration):
    def get_metadata(url, headers=None):
        return MockResponse(WORKFLOW_METADATA, 200)

    monkeypatch.setattr(requests, "get", get_metadata)
    runner = click.testing.CliRunner()
    result = runner.invoke(cli.main, ["metadata", "36"])
    assert "workflowName" in result.output

    def get_tasks(url, headers=None):
        return MockResponse(WORKFLOW_METADATA, 201)


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_cli_submit(configuration, monkeypatch, sample_workflow):
    root = sample_workflow

    wdl = os.path.join(root, "workflow", "sample.wdl")
    inputs = os.path.join(root, "workflow", "sample.json")

    def mock_get(url, headers=None):
        if "user" in url:
            result = {"email": "joe@lbl.gov", "uid": "jdoe", "name": "John Doe"}
        else:
            result = {
                "site_id": "CORI",
                "globus_endpoint": "abcdeqerawr13423sdasd",
                "globus_host_path": "/",
                "uploads_dir": "/global/cscratch1/sd/jaws_jtm/jaws-dev/uploads",
                "max_ram_gb": 1024,
            }
        return MockResponse(result, 200)

    def mock_post(url, data=None, files=None, headers={}):
        return MockResponse(
            {"run_id": "36", "site_id": "CORI", "tag": None, "output_dir": "/a/b/c"},
            201,
        )

    monkeypatch.setattr(requests, "get", mock_get)
    monkeypatch.setattr(requests, "post", mock_post)

    runner = click.testing.CliRunner()
    result = runner.invoke(cli.main, ["submit", wdl, inputs, "CORI"])
    assert result.exit_code == 0


def test_get(configuration, monkeypatch):
    def mock__run_status(run_id, verbose):
        if run_id == "1":
            return {
                "status": "download complete",
                "output_dir": "/data/repo/dir/mockuser/run1",
            }
        else:
            return {"status": "submitted", "output_dir": "/data/repo/dir/mockuser/run1"}

    def mock_run(args, **kwargs):
        if args[0] == "rsync" and args[1] == "-a" and len(args) == 4:
            return
        else:
            raise ValueError

    def mock_rsync(src, dest, options):
        class Result:
            def __init__(self):
                self.returncode = 0

        return Result()

    def mock_get_outputs(run_id, src, dest):
        pass

    def mock_get_complete(run_id, src, dest):
        pass

    def mock_makedirs(path):
        pass
    
    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(cli, "_run_status", mock__run_status)
    monkeypatch.setattr(cli, "_get_outputs", mock_get_outputs)
    monkeypatch.setattr(cli, "_get_complete", mock_get_complete)
    monkeypatch.setattr(os, "makedirs", mock_makedirs)

    from jaws_client import workflow

    monkeypatch.setattr(workflow, "rsync", mock_rsync)

    runner = click.testing.CliRunner()

    # a completed run
    result = runner.invoke(cli.main, ["get", "1", "/home/mockuser/mydir"])
    assert result.exit_code == 0

    # an incomplete run
    result = runner.invoke(cli.main, ["get", "2", "/home/mockuser/mydir"])
    assert result.exit_code != 0


def test_cancel_OK(monkeypatch, configuration):
    """Check if cancel run producs expected JSON output in case of a successful cancel."""

    def put_cancel(url, data={}, files={}, headers={}):
        return MockResponse({"cancel": "OK"}, 200)

    monkeypatch.setattr(requests, "put", put_cancel)
    runner = click.testing.CliRunner()
    result = runner.invoke(cli.main, ["cancel", "35"])
    assert result.exit_code == 0


def test_cancel_ERR(monkeypatch, configuration):
    """Check if cancel run producs expected JSON output in case of an error in cancel run."""

    def put_cancel(url, headers={}):
        err = {
            "error": "That Run had already been cancelled",
            "status": 400,
            "title": "Bad Request",
            "type": "about:blank",
        }
        return MockResponse(err, 400)

    monkeypatch.setattr(requests, "put", put_cancel)
    runner = click.testing.CliRunner()
    result = runner.invoke(cli.main, ["cancel", "35"])
    assert result.exit_code == 1

def test_utc_to_local(monkeypatch):
    from pytz import timezone
    local_time = cli._utc_to_local('2021-07-06 22:42:55')
    assert '2021-07-06 15:42:55' in local_time
