import json
import os
from deepdiff import DeepDiff
from jaws_site import cromwell
from tests.conftest import mock_task_status_table, mock_task_summary_table, this_date

tests_dir = os.path.dirname(os.path.abspath(__file__))

example_cromwell_url = "http://localhost:8000"
example_workflows_url = f"{example_cromwell_url}/api/workflows/v1"
crom = cromwell.Cromwell(example_cromwell_url)

# simple successful run
example_cromwell_run_id_1 = "ee30d68f-39d4-4fde-85c2-afdecce2bad3"

# successful run with subs
example_cromwell_run_id_2 = "c720836c-0931-4ddc-8366-774160e05531"

# simple failed run
example_cromwell_run_id_3 = "dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15"

# failed, missing infile
example_cromwell_run_id_4 = "469bcdd6-d67f-455f-9475-438349f41631"

# successful run with subs
example_cromwell_run_id_5 = "5c7ca4aa-6bc2-4f90-a80c-6232ee64e1f0"

# failed simple run
example_cromwell_run_id_6 = "cbbbc75f-8920-495c-a290-0a1a5f0d1c20"

# failed run, with subs (sub's docker container not found)
example_cromwell_run_id_7 = "dcb85c55-bf74-4f63-bce3-fe61f7b84ebb"

# successful run, with outputs list
example_cromwell_run_id_8 = "5a2cbafe-56ed-42aa-955d-fef8cb5014bb"

# successful AWS run
example_cromwell_run_id_9 = "f4f5afd1-79f5-497a-9612-baed76dc365d"


def __load_example_output_from_file(cromwell_run_id, output_type):
    with open(f"{tests_dir}/{cromwell_run_id}.{output_type}.json", "r") as fh:
        output = json.load(fh)
    return output


def test_get_metadata(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expectedWorkflowRoot = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # noqa
    metadata1 = crom.get_metadata(example_cromwell_run_id_1)
    assert expectedWorkflowRoot == metadata1.get("workflowRoot")


def test_workflow_root(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_9}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_9, "metadata"),
    )
    expectedWorkflowRoot = "s3://jaws-site/cromwell-execution/jgi_meta/f4f5afd1-79f5-497a-9612-baed76dc365d"  # noqa
    metadata = crom.get_metadata(example_cromwell_run_id_9)
    assert expectedWorkflowRoot == metadata.workflow_root()


def test_task(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    task = metadata.tasks["main_workflow.hello_and_goodbye_1"]
    shard_index = -1  # not a sharded task
    attempt = 1  # generally only 1 attempt
    assert task.calls[shard_index][attempt].execution_status == "Done"


def test_task_subworkflow(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata2 = crom.get_metadata(example_cromwell_run_id_2)
    task = metadata2.tasks["main_workflow.hello_and_goodbye_1"]
    shard_index = -1  # not a sharded task
    attempt = 1  # generally only 1 attempt
    taskMetadata = task.calls[shard_index][attempt].subworkflow
    assert taskMetadata.get("rootWorkflowId") == "c720836c-0931-4ddc-8366-774160e05531"


def test_task_summary(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    actual = metadata.task_summary()
    expected = [
        [
            "main_workflow.goodbye",
            "12129",
            False,
            "00:10:00",
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/c720836c-0931-4ddc-8366-774160e05531/call-goodbye/execution",  # noqa
        ],
        [
            "main_workflow.hello",
            "12130",
            False,
            "00:10:00",
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/c720836c-0931-4ddc-8366-774160e05531/call-hello/execution",  # noqa
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "12134",
            False,
            "00:10:00",
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/c720836c-0931-4ddc-8366-774160e05531/call-hello_and_goodbye_1/sub.hello_and_goodbye/6870a657-27df-4972-9465-88d769b81e49/call-goodbye/execution",  # noqa
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "12133",
            False,
            "00:10:00",
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/c720836c-0931-4ddc-8366-774160e05531/call-hello_and_goodbye_1/sub.hello_and_goodbye/6870a657-27df-4972-9465-88d769b81e49/call-hello/execution",  # noqa
        ],
        [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
            "12131",
            False,
            "00:10:00",
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/c720836c-0931-4ddc-8366-774160e05531/call-hello_and_goodbye_2/sub.hello_and_goodbye/5689d65d-51bf-4d7f-b134-cd086ba6195b/call-goodbye/execution",  # noqa
        ],
        [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
            "12132",
            False,
            "00:10:00",
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/c720836c-0931-4ddc-8366-774160e05531/call-hello_and_goodbye_2/sub.hello_and_goodbye/5689d65d-51bf-4d7f-b134-cd086ba6195b/call-hello/execution",  # noqa
        ],
    ]
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_5}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_5, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_5)
    actual = metadata.task_summary()
    expected = __load_example_output_from_file(
        example_cromwell_run_id_5, "task-summary"
    )
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


def test_task_stdout(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_1)
    task = metadata.tasks["fq_count.count_seqs"]
    shard_index = -1  # not a sharded task
    attempt = 1  # generally only 1 attempt
    call = task.calls[shard_index][attempt]
    expected_stderr = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr"  # noqa
    expected_stdout = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout"  # noqa
    assert call.stderr == expected_stderr
    assert call.stdout == expected_stdout


def test_metadata_tasks(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )

    metadata = crom.get_metadata(example_cromwell_run_id_2)
    for task_name, task in metadata.tasks.items():
        assert len(task.data) > 0
        call = task.data[0]
        shard_index = call["shardIndex"]
        attempt = call["attempt"]
        if "subWorkflowMetadata" in call:
            subworkflow = task.subworkflows[shard_index][attempt]
            assert isinstance(subworkflow, cromwell.Metadata)
        else:
            assert "jobId" in call


def test_errors(requests_mock, monkeypatch):
    def mock_read_file(path):
        return None

    monkeypatch.setattr(cromwell, "_read_file", mock_read_file)

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_errors_report_1 = {}
    metadata_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_errors_report_1 = metadata_1.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_1, expected_errors_report_1, ignore_order=True
            )
        )
        is False
    )

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_3}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_3, "metadata"),
    )
    expected_errors_report_3 = __load_example_output_from_file(
        example_cromwell_run_id_3, "errors"
    )
    metadata_3 = crom.get_metadata(example_cromwell_run_id_3)
    actual_errors_report_3 = metadata_3.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_3, expected_errors_report_3, ignore_order=True
            )
        )
        is False
    )

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_4}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_4, "metadata"),
    )
    expected_errors_report_4 = __load_example_output_from_file(
        example_cromwell_run_id_4, "errors"
    )
    metadata_4 = crom.get_metadata(example_cromwell_run_id_4)
    actual_errors_report_4 = metadata_4.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_4, expected_errors_report_4, ignore_order=True
            )
        )
        is False
    )

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_6}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_6, "metadata"),
    )
    expected_errors_report_6 = __load_example_output_from_file(
        example_cromwell_run_id_6, "errors"
    )
    metadata_6 = crom.get_metadata(example_cromwell_run_id_6)
    actual_errors_report_6 = metadata_6.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_6, expected_errors_report_6, ignore_order=True
            )
        )
        is False
    )

    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_7}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_7, "metadata"),
    )
    expected_errors_report_7 = __load_example_output_from_file(
        example_cromwell_run_id_7, "errors"
    )
    metadata_7 = crom.get_metadata(example_cromwell_run_id_7)
    actual_errors_report_7 = metadata_7.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_7, expected_errors_report_7, ignore_order=True
            )
        )
        is False
    )


def test_get_outputs(requests_mock):
    # test 1 : outputs scalar
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_outputs_1 = {
        "fq_count.outfile": "./call-count_seqs/execution/num_seqs.txt"
    }
    ex_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_outputs_1 = ex_1.outputs(relpath=True)
    assert (
        bool(DeepDiff(actual_outputs_1, expected_outputs_1, ignore_order=True)) is False
    )

    # test 2 : outputs list
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_8}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_8, "metadata"),
    )
    expected_outputs_8 = {
        "create_lastdb.dbFiles": [
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.bck",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.des",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.prj",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.sds",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.ssp",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.suf",
            "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.tis",
        ]
    }
    ex_8 = crom.get_metadata(example_cromwell_run_id_8)
    actual_outputs_8 = ex_8.outputs(relpath=True)
    assert (
        bool(DeepDiff(actual_outputs_8, expected_outputs_8, ignore_order=True)) is False
    )


def test_outfiles(requests_mock):
    # test 1 : outputs scalar
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    expected_outfiles_1 = ["./call-count_seqs/execution/num_seqs.txt"]
    ex_1 = crom.get_metadata(example_cromwell_run_id_1)
    actual_outfiles_1 = ex_1.outfiles(relpath=True)
    assert (
        bool(DeepDiff(actual_outfiles_1, expected_outfiles_1, ignore_order=True))
        is False
    )

    # test 2 : outputs list
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_8}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_8, "metadata"),
    )
    expected_outfiles_8 = [
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.bck",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.des",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.prj",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.sds",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.ssp",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.suf",
        "./call-lastdb/execution/glob-457b18ddcc95e1a8de02cd6f6cc84b25/refGenomes.faa.tis",
    ]
    ex_8 = crom.get_metadata(example_cromwell_run_id_8)
    actual_outfiles_8 = ex_8.outfiles(relpath=True)
    assert (
        bool(DeepDiff(actual_outfiles_8, expected_outfiles_8, ignore_order=True))
        is False
    )

    # test 9 : aws outputs
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_9}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_9, "metadata"),
    )
    expected_outfiles_9 = [
        "./call-create_agp/cacheCopy/assembly.contigs.fasta",
        "./call-read_mapping_pairs/covstats.txt",
        "./call-read_mapping_pairs/pairedMapped.sam.gz",
        "./call-bbcms/cacheCopy/resources.log",
        "./call-create_agp/cacheCopy/assembly.scaffolds.fasta",
        "./call-bbcms/cacheCopy/input.corr.fastq.gz",
        "./call-assy/resources.log",
        "./call-read_mapping_pairs/resources.log",
        "./call-assy/spades3/spades.log",
        "./call-create_agp/cacheCopy/resources.log",
        "./call-read_mapping_pairs/pairedMapped_sorted.bam.bai",
        "./call-bbcms/cacheCopy/counts.metadata.json",
        "./call-bbcms/cacheCopy/stdout.log",
        "./call-assy/spades3/scaffolds.fasta",
        "./call-bbcms/cacheCopy/input.corr.left.fastq.gz",
        "./call-bbcms/cacheCopy/stderr.log",
        "./call-bbcms/cacheCopy/readlen.txt",
        "./call-bbcms/cacheCopy/unique31mer.txt",
        "./call-create_agp/cacheCopy/assembly.scaffolds.legend",
        "./call-bbcms/cacheCopy/input.corr.right.fastq.gz",
        "./call-read_mapping_pairs/pairedMapped_sorted.bam",
        "./call-create_agp/cacheCopy/assembly.agp",
    ]
    ex_9 = crom.get_metadata(example_cromwell_run_id_9)
    actual_outfiles_9 = ex_9.outfiles(relpath=True)
    assert (
        bool(DeepDiff(actual_outfiles_9, expected_outfiles_9, ignore_order=True))
        is False
    )


def test_job_summary(monkeypatch):
    def mock_task_summary(self):
        example_task_summary = [
            ["main_workflow.goodbye", "12129", False, "0:10:00", "Success", "/a/b/c"],
            ["main_workflow.hello", "12130", False, "0:10:00", "Success", "/a/b/c"],
            [
                "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
                "12134",
                False,
                "0:10:00",
                "Success",
                "/a/b/c",
            ],
            [
                "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
                "12133",
                False,
                "0:10:00",
                "Success",
                "/a/b/c",
            ],
            [
                "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
                "12131",
                False,
                "0:10:00",
                "Success",
                "/a/b/c",
            ],
            [
                "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
                "12132",
                False,
                "0:10:00",
                "Success",
                "/a/b/c",
            ],
        ]
        return example_task_summary

    monkeypatch.setattr(cromwell.Metadata, "task_summary", mock_task_summary)

    expected_result = {
        "12129": ["main_workflow.goodbye", "0:10:00"],
        "12130": ["main_workflow.hello", "0:10:00"],
        "12134": [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "0:10:00",
        ],
        "12133": [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "0:10:00",
        ],
        "12131": [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
            "0:10:00",
        ],
        "12132": [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
            "0:10:00",
        ],
    }

    mock_data = {"id": "EXAMPLE_CROMWELL_RUN_ID"}
    metadata = cromwell.Metadata(mock_data)
    actual_result = metadata.job_summary()

    assert bool(DeepDiff(actual_result, expected_result, ignore_order=True)) is False

def test_get_errors(monkeypatch):

    example_cromwell_run_id = "AAAA-BBBB-CCCC"
    example_cromwell_errors_report = {
        "calls": {
            "example_sub": [
                {
                    "subWorkflowMetadata": {
                        "calls": {
                            "example_task": [
                                {
                                    "jobId": "300",
                                }
                            ],
                        }
                    }
                }
            ],
        }
    }
    example_task_log_errors = {
        "300": ["Out of time!"],
    }
    expected_errors_report = {
        "calls": {
            "example_sub": [
                {
                    "subWorkflowMetadata": {
                        "calls": {
                            "example_task": [
                                {
                                    "jobId": "300",
                                    "taskLog": ["Out of time!"],
                                }
                            ],
                        }
                    }
                }
            ],
        }
    }

    def mock_get_cromwell_errors_report(cromwell_run_id):
        return example_cromwell_errors_report

    monkeypatch.setattr(
        errors, "_get_cromwell_errors_report", mock_get_cromwell_errors_report
    )

    def mock_select_task_log_error_messages(session, cromwell_job_id):
        if cromwell_job_id in example_task_log_errors:
            return example_task_log_errors[cromwell_job_id]
        else:
            return []

    monkeypatch.setattr(
        tasks, "_select_task_log_error_messages", mock_select_task_log_error_messages
    )

    mock_session = None
    actual_errors_report = errors.get_errors(mock_session, example_cromwell_run_id)
    assert (
        bool(DeepDiff(actual_errors_report, expected_errors_report, ignore_order=True))
        is False
    )

def test_task_status_table(monkeypatch):
    def mock_task_log(self):
        example_log = [
            [
                "main.ex_task_1",
                "2222",
                False,
                "created",
                "ready",
                "2020-03-22 12:47:10",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                False,
                "ready",
                "queued",
                "2020-03-22 12:47:20",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                False,
                "queued",
                "pending",
                "2020-03-22 12:47:25",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                False,
                "pending",
                "running",
                "2020-03-22 12:48:01",
                None,
            ],
            [
                "main.ex_task_2",
                "2223",
                False,
                "created",
                "ready",
                "2020-03-22 12:48:11",
                None,
            ],
            [
                "main.ex_task_2",
                "2223",
                False,
                "ready",
                "queued",
                "2020-03-22 12:48:16",
                None,
            ],
        ]
        return example_log

    example_run_id = 1
    expected_status = [
        [
            "main.ex_task_1",
            "2222",
            False,
            "running",
            "2020-03-22 12:48:01",
            None,
        ],
        [
            "main.ex_task_2",
            "2223",
            False,
            "queued",
            "2020-03-22 12:48:16",
            None,
        ],
    ]

    monkeypatch.setattr(TaskLog, "task_log", mock_task_log)
    mock_session = None

    tasks = TaskLog(mock_session, run_id=example_run_id)
    task_status = tasks.task_status_table()
    assert bool(DeepDiff(task_status, expected_status, ignore_order=True)) is False


def test_run_status(monkeypatch):
    def mock_task_status(session):
        example_status = []
        if run_id == 102:
            example_status = [
                [
                    "main.ex_task_1",
                    "2222",
                    False,
                    "ready",
                    "2020-03-22 12:47:10",
                    None,
                ],
            ]
        elif run_id == 103:
            example_status = [
                [
                    "main.ex_task_2",
                    "2222",
                    False,
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ]
            ]
        elif run_id == 104:
            example_status = [
                [
                    "main.ex_task_2",
                    "2222",
                    False,
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ],
                [
                    "main.ex_task_3",
                    "2222",
                    False,
                    "pending",
                    "2020-03-22 12:47:25",
                    None,
                ],
            ]
        elif run_id == 105:
            example_status = [
                [
                    "main.ex_task_4",
                    "2222",
                    False,
                    "running",
                    "2020-03-22 12:48:01",
                    None,
                ],
                [
                    "main.ex_task_5",
                    "2223",
                    False,
                    "ready",
                    "2020-03-22 12:48:11",
                    None,
                ],
            ]
        elif run_id == 106:
            example_status = [
                [
                    "main.ex_task_6",
                    "2223",
                    False,
                    "success",
                    "2020-03-22 12:48:16",
                    None,
                ],
            ]
        return example_status

    run_id_and_expected = {
        101: None,
        102: "queued",
        103: "queued",
        104: "queued",
        105: "running",
        106: "running",
    }

    monkeypatch.setattr(TaskLog, "task_status", mock_task_status)
    mock_session = None

    for run_id, expected in run_id_and_expected.items():
        assert tasks.get_run_status(mock_session, run_id) == expected


def test_task_summary_table(monkeypatch):
    def mock_task_log(self):
        self._task_log = [
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "created",
                "ready",
                "2021-12-07 20:39:09",
                None,
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "ready",
                "queued",
                "2021-12-07 20:39:09",
                None,
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "queued",
                "pending",
                "2021-12-07 20:39:09",
                "slurm_jid=45352308",
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "pending",
                "running",
                "2021-12-07 20:39:16",
                None,
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "running",
                "success",
                "2021-12-07 20:39:16",
                None,
            ],
        ]
        return self._task_log

    monkeypatch.setattr(TaskLog, "task_log", mock_task_log)

    def mock_get_cromwell_run_id(self):
        self._cromwell_run_id = "AAAA"
        return self._cromwell_run_id

    monkeypatch.setattr(TaskLog, "_get_cromwell_run_id", mock_get_cromwell_run_id)

    def mock_cromwell_job_summary(self):
        self._cromwell_job_summary = {"8919": ["fq_count.count_seqs", "00:10:00"]}
        return self._cromwell_job_summary

    monkeypatch.setattr(TaskLog, "cromwell_job_summary", mock_cromwell_job_summary)

    expected = [
        [
            "fq_count.count_seqs",
            "8919",
            False,
            "success",
            "2021-12-07 20:39:09",
            "0:00:07",
            "0:00:00",
            "00:10:00",
        ]
    ]

    example_run_id = 1
    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    actual = tasks.task_summary_table()
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False


def test_get_task_cromwell_dir_mapping(monkeypatch):
    class TaskMetadata:
        def __init__(self, task_data):
            self.tasks = task_data

        def summary(self):
            return self.tasks

    def mock_cromwell(*args, **kwargs):
        class CromMetadata:
            task1 = [
                [
                    "align.stats",
                    1,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-stats/execution",
                ],
            ]
            task2 = [
                [
                    "align.shard_wf:shard_wf.indexing",
                    2,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-indexing/execution",
                ],
                [
                    "align.shard_wf:shard_wf.map[0]",
                    3,
                    False,
                    "01:00:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-map/shard-0/execution",
                ],
                [
                    "align.shard_wf:shard_wf.merge",
                    4,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-merge/execution",
                ],
                [
                    "align.shard_wf:shard_wf.shard",
                    5,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-shard/execution",
                ],
            ]
            task3 = [
                [
                    "no_cromwell_dir",
                    1,
                    False,
                    "00:30:00",
                    "",
                ],
            ]

            def __init__(self):
                self.tasks = {
                    "align.stats": TaskMetadata(CromMetadata.task1),
                    "align.bbmap_shard_wf": TaskMetadata(CromMetadata.task2),
                    "no_cromwell_dir": TaskMetadata(CromMetadata.task3),
                }

        return CromMetadata()

    monkeypatch.setattr(cromwell.Cromwell, "get_metadata", mock_cromwell)

    expected = {
        "cromwell-executions/align/C1/call-stats/execution": "align.stats",
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-indexing/execution": "align.shard_wf:shard_wf.indexing",  # noqa
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-map/shard-0/execution": "align.shard_wf:shard_wf.map[0]",  # noqa
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-merge/execution": "align.shard_wf:shard_wf.merge",  # noqa
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-shard/execution": "align.shard_wf:shard_wf.shard",  # noqa
    }
    mock_session = None
    tasks = TaskLog(mock_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    actual = tasks.get_task_cromwell_dir_mapping()
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False


def test_task_summary(mock_db_session, monkeypatch):
    monkeypatch.setattr(TaskLog, "task_summary_table", mock_task_summary_table)

    exp_results = {
        "task_abcd": {
            "cached": False,
            "cromwell_job_id": "cromwell_abcd",
            "max_time": "03:00:00",
            "queue_wait": "01:00:00",
            "queued": this_date,
            "result": "success",
            "run_time": "02:00:00",
        },
        "task_efgh": {
            "cached": True,
            "cromwell_job_id": "cromwell_efgh",
            "max_time": "06:00:00",
            "queue_wait": "04:00:00",
            "queued": this_date,
            "result": "success",
            "run_time": "05:00:00",
        },
    }

    tasks = TaskLog(mock_db_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    obs_results = tasks.task_summary()
    assert bool(DeepDiff(obs_results, exp_results, ignore_order=True)) is False


def test_task_status(mock_db_session, monkeypatch):
    monkeypatch.setattr(TaskLog, "task_status_table", mock_task_status_table)

    exp_results = {
        "task_abcd": {
            "cromwell_job_id": "cromwell_abcd",
            "reason": "reason_abcd",
            "status": "success",
            "timestamp": this_date,
            "cached": False,
        },
        "task_efgh": {
            "cromwell_job_id": "cromwell_efgh",
            "reason": "reason_efgh",
            "status": "success",
            "timestamp": this_date,
            "cached": True,
        },
    }

    tasks = TaskLog(mock_db_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    obs_results = tasks.task_status()
    assert bool(DeepDiff(obs_results, exp_results, ignore_order=True)) is False
