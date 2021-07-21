import pytest
import os
import shutil
from jaws_client import api
from jaws_client.api import Client


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
    "runblastplus_sub.task1\t1\t43\tready\tqueued\t2020-06-10 13:42:44\t\tThe job was received by JTM-manager and sent to JTM-worker\n"  # noqa
    "runblastplus_sub.task2\t1\t44\tqueued\tpending\t2020-06-10 13:43:36\t\tThe job was receive by JTM-worker and is awaiting resources\n"  # noqa
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


QUEUE = [
    [
        35,
        "2020-04-04T16:55:35Z",
        "running",
        "65f2f4df-2a6c-4881-a3b0-3141107ac668",
        "1aba54ac-7695-11ea-9615-0afc9e7dd773",
    ]
]

TASK_STATUS_JSON = [
    ["bbtools.alignment", 1, 432, "Queued", "Running", "2020-04-03 20:32:52", ""]
]

TASK_STATUS_TEXT = (
    "#TASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON\n"
    "bbtools.alignment\t1\t432\tqueued\trunning\t2020-04-03 20:32:52\t\n"
)


def test_queue(monkeypatch, configuration):
    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        result = [
            [
                "36",
                "2020-04-04T17:56:52Z",
                "running",
                "bafa20aa-59b3-4d33-9f4c-fb30696e324f",
                "aa791756-769d-11ea-af53-0201714f6eab",
            ],
        ]
        return result

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()
    result = jaws.queue()
    assert result[0][2] == "running"


def test_history(monkeypatch, configuration):
    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        return HISTORY

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()
    result = jaws.history()

    expected = {0: "33", 1: "34", 2: "35", 3: "36"}

    for index, run_id in expected.items():
        assert result[index][0] == run_id


def test_status(monkeypatch, configuration):
    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        job_status = {"status": "Running"}
        if url.endswith("/complete"):
            job_status["output_dir"] = "/foo/bar"
        return job_status

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()

    result = jaws.status("36", True)
    assert result["status"] == "Running"
    assert result["output_dir"] == "/foo/bar"

    result = jaws.status("36")
    assert result["status"] == "Running"
    assert "output_dir" not in result


def test_metadata(monkeypatch, configuration):
    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        return WORKFLOW_METADATA

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()
    result = jaws.metadata("36")
    assert "workflowName" in result


@pytest.mark.skipif(
    shutil.which("womtool") is None, reason="WOMTool needs to be installed."
)
def test_submit(configuration, monkeypatch, sample_workflow):
    root = sample_workflow

    wdl = os.path.join(root, "workflow", "sample.wdl")
    inputs = os.path.join(root, "workflow", "sample.json")

    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        if rest_op == "GET":
            if "user" in url:
                return {"email": "joe@lbl.gov", "uid": "jdoe", "name": "John Doe"}
            else:
                return {
                    "site_id": "CORI",
                    "globus_endpoint": "abcdeqerawr13423sdasd",
                    "globus_host_path": "/",
                    "uploads_dir": "/global/cscratch1/sd/jaws_jtm/jaws-dev/uploads",
                    "max_ram_gb": 1024,
                }
        else:
            return {
                "run_id": "36",
                "site_id": "CORI",
                "tag": None,
                "output_dir": "/a/b/c",
            }

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()
    result = jaws.submit(wdl, inputs, "CORI")
    assert "run_id" in result


def test_get(configuration, monkeypatch):
    def mock_status(self, run_id: int, verbose: bool = False) -> dict:
        if run_id == "1":
            return {
                "status": "download complete",
                "output_dir": "/data/repo/dir/mockuser/run1",
            }
        else:
            return {"status": "submitted", "output_dir": "/data/repo/dir/mockuser/run1"}

    monkeypatch.setattr(api.Client, "status", mock_status)

    def mock_rsync(src, dest, options):
        class Result:
            def __init__(self):
                self.returncode = 0

        return Result()

    from jaws_client import workflow

    monkeypatch.setattr(workflow, "rsync", mock_rsync)

    jaws = Client()

    # a completed run
    jaws.get("1", "/home/mockuser/mydir")

    # an incomplete run
    with pytest.raises(api.JawsClientError):
        jaws.get("2", "/home/mockuser/mydir")


def test_cancel_OK(monkeypatch, configuration):
    """Check if cancel run producs expected JSON output in case of a successful cancel."""

    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        return {"cancel": "OK"}

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()
    result = jaws.cancel("35")
    assert result["cancel"] == "OK"


def test_cancel_ERR(monkeypatch, configuration):
    """Check if cancel run producs expected JSON output in case of an error in cancel run."""

    def mock_request(self, rest_op, url, data={}, files={}) -> dict:
        return {
            "error": "That Run had already been cancelled",
            "status": 400,
            "title": "Bad Request",
            "type": "about:blank",
        }

    monkeypatch.setattr(api.Client, "_request", mock_request)

    jaws = Client()
    result = jaws.cancel("35")
    assert "error" in result
