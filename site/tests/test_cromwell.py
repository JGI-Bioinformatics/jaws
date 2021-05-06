import json
import os
from deepdiff import DeepDiff
from jaws_site import cromwell

tests_dir = os.path.dirname(os.path.abspath(__file__))

example_cromwell_url = "http://localhost:8000"
example_workflows_url = f"{example_cromwell_url}/api/workflows/v1"
crom = cromwell.Cromwell(example_cromwell_url)

example_cromwell_run_id_1 = (
    "ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # simple successful run
)
example_cromwell_run_id_2_main = (
    "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a"  # successful run, with subs
)
example_cromwell_run_id_2_sub_1 = "7408a4f1-bc85-49ba-8d5f-c886261ab6a0"
example_cromwell_run_id_2_sub_2 = "89d86efc-dd04-48aa-a65f-21fb9d0c8be3"
example_cromwell_run_id_3 = "dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15"  # simple failed run
example_cromwell_run_id_4 = (
    "469bcdd6-d67f-455f-9475-438349f41631"  # failed, missing infile
)
example_cromwell_run_id_5 = (
    "36f86fab-830c-4022-86be-a175c61989d7"  # failed run, with scatter
)


def __load_example_metadata_from_file(cromwell_run_id):
    with open(f"{tests_dir}/{cromwell_run_id}.json", "r") as fh:
        metadata = json.load(fh)
    return metadata


def test_metadata(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_1),
    )
    expectedWorkflowRoot = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # noqa

    metadata = crom.get_metadata(example_cromwell_run_id_1)
    workflowRoot = metadata.get("workflowRoot")
    assert workflowRoot == expectedWorkflowRoot


def test_task(requests_mock):
    example_metadata = __load_example_metadata_from_file(example_cromwell_run_id_1)
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=example_metadata,
    )
    example_task_name = "fq_count.count_seqs"
    example_calls = example_metadata["calls"][example_task_name]

    task = cromwell.Task(example_workflows_url, example_task_name, example_calls)
    assert task.name == example_task_name
    assert task.is_subworkflow() is False
    assert task.calls[0]["attempt"] == 1
    assert task.calls[0]["jobId"] == "30"
    assert (
        task.get("callRoot", 1)
        == "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs"  # noqa
    )  # noqa


def test_task_subworkflow(requests_mock):
    example_metadata = __load_example_metadata_from_file(example_cromwell_run_id_2_main)
    example_metadata_sub_1 = __load_example_metadata_from_file(
        example_cromwell_run_id_2_sub_1
    )
    example_metadata_sub_2 = __load_example_metadata_from_file(
        example_cromwell_run_id_2_sub_2
    )
    assert example_metadata
    assert example_metadata_sub_1
    assert example_metadata_sub_2
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_main}/metadata",
        json=example_metadata,
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_1}/metadata",
        json=example_metadata_sub_1,
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_2}/metadata",
        json=example_metadata_sub_2,
    )

    example_task_name = "main_workflow.hello_and_goodbye_1"  # task in sub
    example_calls = example_metadata["calls"][example_task_name]

    task = cromwell.Task(example_workflows_url, example_task_name, example_calls)
    expected_subworkflow_job_ids = (5482, 5484)
    for call in task.calls:
        if "jobId" in call:
            assert call["jobId"] in expected_subworkflow_job_ids
        else:
            assert call["subWorkflowId"] == example_cromwell_run_id_2_sub_1


def test_get_all_metadata(requests_mock):
    """Given the uuid of a workflow with a subworkflow, we expect a dict with { uuid => metadata }"""
    example_metadata = __load_example_metadata_from_file(example_cromwell_run_id_2_main)
    example_metadata_sub_1 = __load_example_metadata_from_file(
        example_cromwell_run_id_2_sub_1
    )
    example_metadata_sub_2 = __load_example_metadata_from_file(
        example_cromwell_run_id_2_sub_2
    )
    assert example_metadata
    assert example_metadata_sub_1
    assert example_metadata_sub_2
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_main}/metadata",
        json=example_metadata,
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_1}/metadata",
        json=example_metadata_sub_1,
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_2}/metadata",
        json=example_metadata_sub_2,
    )

    c = cromwell.Cromwell(example_cromwell_url)
    result = c.get_all_metadata(example_cromwell_run_id_2_main)
    assert (
        bool(
            DeepDiff(
                result[example_cromwell_run_id_2_main],
                example_metadata,
                ignore_order=True,
            )
        )
        is False
    )
    assert (
        bool(
            DeepDiff(
                result[example_cromwell_run_id_2_sub_1],
                example_metadata_sub_1,
                ignore_order=True,
            )
        )
        is False
    )
    assert (
        bool(
            DeepDiff(
                result[example_cromwell_run_id_2_sub_2],
                example_metadata_sub_2,
                ignore_order=True,
            )
        )
        is False
    )


def test_task_summary(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_main}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_2_main),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_1}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_2_sub_1),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_2}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_2_sub_2),
    )

    metadata = crom.get_metadata(example_cromwell_run_id_2_main)
    result = metadata.task_summary()
    expected = [
        [example_cromwell_run_id_2_main, "main_workflow.goodbye", 1, "5480"],
        [example_cromwell_run_id_2_main, "main_workflow.hello", 1, "5481"],
        [example_cromwell_run_id_2_sub_2, "hello_and_goodbye.goodbye", 1, "5483"],
        [example_cromwell_run_id_2_sub_2, "hello_and_goodbye.hello", 1, "5485"],
        [example_cromwell_run_id_2_sub_1, "hello_and_goodbye.goodbye", 1, "5482"],
        [example_cromwell_run_id_2_sub_1, "hello_and_goodbye.hello", 1, "5484"],
    ]
    assert bool(DeepDiff(result, expected, ignore_order=True)) is False


def test_task_stdout(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_1),
    )

    metadata = crom.get_metadata(example_cromwell_run_id_1)
    expected_stderr = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr"  # noqa
    expected_stdout = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout"  # noqa
    assert len(metadata.tasks) == 1
    task = metadata.tasks[0]
    assert task.stderr() == expected_stderr
    assert task.stdout() == expected_stdout


def test_metadata_tasks(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_main}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_2_main),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_1}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_2_sub_1),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_2_sub_2}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_2_sub_2),
    )

    metadata = crom.get_metadata(example_cromwell_run_id_2_main)
    assert metadata.is_subworkflow() is False
    for task in metadata.tasks:
        assert task.name is not None
        assert len(task.calls) > 0
        for call in task.calls:
            assert "attempt" in call
            assert call["attempt"] == 1
            assert "jobId" in call or "subWorkflowId" in call
            if "subWorkflowId" in call:
                assert call["subWorkflowId"] in task.subworkflows


def test_errors(requests_mock):
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_1}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_1),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_3}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_3),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_4}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_4),
    )
    requests_mock.get(
        f"{example_cromwell_url}/api/workflows/v1/{example_cromwell_run_id_5}/metadata",
        json=__load_example_metadata_from_file(example_cromwell_run_id_5),
    )

    expected_errors_report_1 = {}

    expected_errors_report_3 = {
        "calls": {
            "fq_count.count_seqs": [
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Unable to start job. Check the stderr file for possible errors: /global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/execution/stderr.submit",  # noqa
                        }
                    ],
                    "jobId": "9999",
                    "runtimeAttributes": {
                        "account": "fungalp",
                        "cluster": "cori",
                        "constraint": "haswell",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "genepool_special",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "stderr": None,
                    "stderr.submit": None,
                }
            ]
        }
    }

    expected_errors_report_4 = {
        "failures": [
            {
                "causedBy": [
                    {
                        "causedBy": [],
                        "message": "Required workflow input 'fq_count.fastq_file' not specified",
                    }
                ],
                "message": "Workflow input processing failed",
            }
        ],
        "inputs": {},
    }

    expected_errors_report_5 = {
        "calls": {
            "jgi_dap_leo.assignGenes": [
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.assignGenes:1:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69550",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 1,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.assignGenes:70:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69454",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 70,
                    "stderr": None,
                    "stderr.submit": None,
                },
            ],
            "jgi_dap_leo.copyOutput_expt": [
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.copyOutput_expt:11:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69694",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 11,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.copyOutput_expt:15:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69670",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 15,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.copyOutput_expt:21:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69646",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 21,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.copyOutput_expt:84:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69502",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 84,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.copyOutput_expt:88:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69622",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 88,
                    "stderr": None,
                    "stderr.submit": None,
                },
            ],
            "jgi_dap_leo.dapStats": [
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.dapStats:50:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69526",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 50,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.dapStats:52:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69574",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 52,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.dapStats:55:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69598",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 55,
                    "stderr": None,
                    "stderr.submit": None,
                },
            ],
            "jgi_dap_leo.findPeaks": [
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.findPeaks:59:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69382",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 59,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.findPeaks:83:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69406",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 83,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.findPeaks:87:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69478",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 87,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.findPeaks:89:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69358",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 89,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.findPeaks:91:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69430",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 91,
                    "stderr": None,
                    "stderr.submit": None,
                },
            ],
            "jgi_dap_leo.trimAlign_expt": [
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.trimAlign_expt:4:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69238",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 4,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.trimAlign_expt:24:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69310",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 24,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.trimAlign_expt:61:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69262",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 61,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.trimAlign_expt:63:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69286",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 63,
                    "stderr": None,
                    "stderr.submit": None,
                },
                {
                    "attempt": 1,
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Job jgi_dap_leo.trimAlign_expt:69:1 exited with return code 79 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.",  # noqa
                        }
                    ],
                    "jobId": "69334",
                    "runtimeAttributes": {
                        "account": "jgi",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "32",
                        "docker": "jfroula/dap_py:3.2",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "memory": "115 GB",
                        "node": "2",
                        "nwpn": "12",
                        "partition": "jgi",
                        "poolname": "dapseq_leo",
                        "qos": "normal",
                        "shared": "0",
                        "time": "05:00:00",
                    },
                    "shardIndex": 69,
                    "stderr": None,
                    "stderr.submit": None,
                },
            ],
        }
    }

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

    metadata_5 = crom.get_metadata(example_cromwell_run_id_5)
    actual_errors_report_5 = metadata_5.errors()
    assert (
        bool(
            DeepDiff(
                actual_errors_report_5, expected_errors_report_5, ignore_order=True
            )
        )
        is False
    )
