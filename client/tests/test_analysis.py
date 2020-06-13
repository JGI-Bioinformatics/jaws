import pytest
import os
import requests
import shutil
import json
from datetime import datetime

import click.testing

from jaws_client.analysis import run
import jaws_client.user


HISTORY = [
    [
        "33",
        "2020-04-03T22:19:45Z",
        "upload failed",
        "c49043b3-4ec7-4874-b09d-c0ca54a74870",
        "null",
    ],
    [
        "34",
        "2020-04-03T22:27:00Z",
        "upload failed",
        "64864a67-04d3-47c9-8cf0-4cb6e0ab020a",
        "null",
    ],
    [
        "35",
        "2020-04-04T16:55:35Z",
        "failed",
        "65f2f4df-2a6c-4881-a3b0-3141107ac668",
        "1aba54ac-7695-11ea-9615-0afc9e7dd773",
    ],
    [
        "36",
        "2020-04-04T17:56:52Z",
        "running",
        "bafa20aa-59b3-4d33-9f4c-fb30696e324f",
        "aa791756-769d-11ea-af53-0201714f6eab",
    ],
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
        "inputs": '{"bbtools.reads":"/global/dna/shared/data/jfroula/JAWS/data/5min_reads.fq","bbtools.ref":"/global/dna/shared/data/jfroula/JAWS/data/5min_ref.fasta"}',  # noqa
        "labels": '{"username": "mamelara"}',
        "options": "{\n\n}",
        "root": "",
        "workflow": 'workflow bbtools { File reads\n    File ref\n\n    call alignment {\n       input: fastq=reads,\n              fasta=ref\n    }\n    call samtools {\n       input: sam=alignment.sam\n   }\n}\n\ntask alignment {\n    File fastq\n    File fasta\n\n    command <<<\n        shifterimg pull jfroula/bbtools:1.2.1 && \\\n        shifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=${fastq} ref=${fasta} out=test.sam\n    >>>\n    output {\n       File sam = "test.sam"\n    }\n}\n\ntask samtools {\n    File sam\n\n    command {\n       shifter --image=jfroula/bbtools:1.2.1 samtools view -b -F0x4 ${sam} | shifter --image=jfroula/bbtools:1.2.1 samtools sort - > test.sorted.bam\n    }\n    output {\n       File bam = "test.sorted.bam"\n    }\n\t#runtime {\n\t#  docker: "jfroula/bbtools:1.2.1"\n\t#}\n}\n',  # noqa
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


RUN_LOG = [
    ["created", "uploading", "2020-06-11 10:01:14", ""],
    ["uploading", "upload complete", "2020-06-11 10:04:23", ""],
]

SUBMISSION = {
    "output_dir": "/global/homes/m/mamelara/out",
    "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    "run_id": 35,
    "site_id": "NERSC",
    "status": "uploading",
    "submission_id": "65f2f4df-2a6c-4881-a3b0-3141107ac668",
    "upload_task_id": "1aba54ac-7695-11ea-9615-0afc9e7dd773",
}


QUEUE = [
    [
        35,
        "2020-04-04T16:55:35Z",
        "running",
        "65f2f4df-2a6c-4881-a3b0-3141107ac668",
        "1aba54ac-7695-11ea-9615-0afc9e7dd773",
    ]
]

TASK_STATUS = [["bbtools.alignment", "Running", "2020-04-03T20:32:52.938Z", ""]]


class MockUser:
    def __init__(self):
        pass

    def header(self):
        return {"Authorization": "Bearer ABCDEFGHIJKLMNOP"}


class MockResult:
    def __init__(self, response, status_code):
        self.response = response
        self.status_code = status_code

    def json(self):
        return self.response

    @property
    def text(self):
        return str(self.response)


@pytest.fixture()
def mock_user(monkeypatch):
    monkeypatch.setattr(jaws_client.user, "User", MockUser)


def test_cli_queue(mock_user, monkeypatch, configuration):
    def get_queue(url, headers=None):
        body = {
            "results": [
                [
                    "36",
                    "2020-04-04T17:56:52Z",
                    "running",
                    "bafa20aa-59b3-4d33-9f4c-fb30696e324f",
                    "aa791756-769d-11ea-af53-0201714f6eab",
                ]
            ]
        }
        return MockResult(body, 200)

    monkeypatch.setattr(requests, "get", get_queue)

    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["queue"])
    assert result.exit_code == 0
    assert "running" in result.output


def test_cli_history(mock_user, monkeypatch, configuration):
    def get_history(url, headers=None):
        body = HISTORY
        return MockResult({"history": body}, 200)

    monkeypatch.setattr(requests, "get", get_history)
    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["history"])
    assert result.exit_code == 0

    for task_id in ["33", "34", "35", "36"]:
        assert task_id in result.output


def test_cli_status(mock_user, monkeypatch, configuration):
    def mock_status_get(url, headers={}):
        """
        Returns a list of attributes from a run ordered as the following:
        run id, submission date, submission id, upload id
        """
        job_status = {"status": "Running"}
        return MockResult(job_status, 200)

    monkeypatch.setattr(requests, "get", mock_status_get)
    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["status", "36"])
    assert result.exit_code == 0
    assert "Running" in result.output


def test_cli_tasks(monkeypatch, mock_user, configuration):
    def get_tasks(url, headers=None):
        return MockResult(TASK_STATUS, 200)

    monkeypatch.setattr(requests, "get", get_tasks)
    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["tasks", "36"])
    assert "Running" in result.output


def test_cli_metadata(monkeypatch, mock_user, configuration):
    def get_metadata(url, headers=None):
        return MockResult(WORKFLOW_METADATA, 200)

    monkeypatch.setattr(requests, "get", get_metadata)
    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["metadata", "36"])
    assert "workflowName" in result.output

    def get_tasks(url, headers=None):
        return MockResult(WORKFLOW_METADATA, 201)


def test_cli_log(monkeypatch, mock_user, configuration):
    def get_log(url, headers=None):
        return MockResult(RUN_LOG, 200)

    valid_states = ['created', 'uploading', 'upload complete',
                    'invalid input', 'submitted', 'submission failed',
                    'queued', 'running', 'succeeded', 'failed', 'cancelling',
                    'cancelled', 'ready', 'downloading', 'download failed',
                    'download complete', 'finished']

    monkeypatch.setattr(requests, "get", get_log)
    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["log", "36"])
    output = result.output.strip()
    output = output.replace("'", '"')  # Convert to valid JSON hack
    run_log = json.loads(output)

    for log_entry in run_log:
        assert log_entry[0] in valid_states
        assert log_entry[1] in valid_states
        assert datetime.strptime(log_entry[2], "%Y-%m-%d %H:%M:%S")


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_cli_submit(configuration, mock_user, monkeypatch, sample_workflow):
    root = sample_workflow

    wdl = os.path.join(root, "workflow", "sample.wdl")
    inputs = os.path.join(root, "workflow", "sample.json")
    subdir = os.path.join(root, "globus", "staging")

    def get_site(url, headers=None):
        body = {
            "globus_basepath": "/NERSC/globus",
            "staging_subdir": "/NERSC/globus/staging",
            "max_ram_gb": "256",
        }
        return MockResult(body, 200)

    def mock_post(url, data=None, files=None, headers={}):
        return MockResult({"run_id": "36"}, 201)

    monkeypatch.setattr(requests, "get", get_site)
    monkeypatch.setattr(requests, "post", mock_post)

    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["submit", wdl, inputs, subdir, "NERSC"])
    assert result.exit_code == 0
