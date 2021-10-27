import json
import os
from deepdiff import DeepDiff
from jaws_site import cromwell

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


def test_task(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata")
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    task = metadata.tasks["main_workflow.hello_and_goodbye_1"]
    assert task.get("executionStatus", -1, 1) == "Done"


def test_task_subworkflow(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata")
    )
    metadata2 = crom.get_metadata(example_cromwell_run_id_2)
    task = metadata2.tasks["main_workflow.hello_and_goodbye_1"]
    taskMetadata = task.get("subWorkflowMetadata")
    assert taskMetadata.get("rootWorkflowId") == "c720836c-0931-4ddc-8366-774160e05531"


def test_task_summary(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_2, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_2)
    actual = metadata.task_summary()
    expected = [
        ["main_workflow.goodbye", "12129", False],
        ["main_workflow.hello", "12130", False],
        ["main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye", "12134", False],
        ["main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello", "12133", False],
        ["main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye", "12131", False],
        ["main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello", "12132", False],
    ]
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_5}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_5, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_5)
    actual = metadata.task_summary()
    print(json.dumps(actual, indent=4, sort_keys=True))  # DEBUG
    expected = __load_example_output_from_file(example_cromwell_run_id_5, "task-summary")
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


def test_task_stdout(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_output_from_file(example_cromwell_run_id_1, "metadata"),
    )
    metadata = crom.get_metadata(example_cromwell_run_id_1)
    task = metadata.tasks["fq_count.count_seqs"]
    expected_stderr = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr"  # noqa
    expected_stdout = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout"  # noqa
    assert task.get("stderr") == expected_stderr
    assert task.get("stdout") == expected_stdout


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
    expected_errors_report_3 = __load_example_output_from_file(example_cromwell_run_id_3, "errors")
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
    expected_errors_report_4 = __load_example_output_from_file(example_cromwell_run_id_4, "errors")
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
    expected_errors_report_6 = __load_example_output_from_file(example_cromwell_run_id_6, "errors")
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
    expected_errors_report_7 = __load_example_output_from_file(example_cromwell_run_id_7, "errors")
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
