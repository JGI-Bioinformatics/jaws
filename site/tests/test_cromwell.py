from deepdiff import DeepDiff
from jaws_site import cromwell

WORKFLOW_ID_EX1 = "ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # simple successful run
WORKFLOW_ID_EX2_MAIN = (
    "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a"  # successful run, with subworkflows
)
WORKFLOW_ID_EX2_SUB1 = "7408a4f1-bc85-49ba-8d5f-c886261ab6a0"
WORKFLOW_ID_EX2_SUB2 = "89d86efc-dd04-48aa-a65f-21fb9d0c8be3"
WORKFLOW_ID_EX3 = "dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15"  # simple failed run

METADATA = {
    "ee30d68f-39d4-4fde-85c2-afdecce2bad3": {  # METADATA FOR WORKFLOW_ID_EX1
        "actualWorkflowLanguage": "WDL",
        "actualWorkflowLanguageVersion": "draft-2",
        "calls": {
            "fq_count.count_seqs": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs",  # noqa
                    "commandLine": 'wc -l /global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/inputs/-230500445/tiny.fastq | perl -ne \'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, " sequences\\n"} else {print STDERR "Invalid Fastq file\\n"}\' > num_seqs.txt',  # noqa
                    "end": "2020-06-10T03:43:54.794Z",
                    "executionEvents": [
                        {
                            "description": "RunningJob",
                            "endTime": "2020-06-10T03:43:54.673Z",
                            "startTime": "2020-06-10T03:38:19.896Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-06-10T03:38:19.896Z",
                            "startTime": "2020-06-10T03:38:19.887Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-06-10T03:38:19.887Z",
                            "startTime": "2020-06-10T03:38:19.887Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-06-10T03:38:19.887Z",
                            "startTime": "2020-06-10T03:38:19.089Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-06-10T03:43:54.794Z",
                            "startTime": "2020-06-10T03:43:54.673Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-06-10T03:38:19.089Z",
                            "startTime": "2020-06-10T03:38:19.089Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {
                        "infile": "/global/cscratch1/sd/jaws/test/uploads/ekirton/NERSC/global/cfs/projectdirs/jaws/test/tiny.fastq"  # noqa
                    },
                    "jobId": "30",
                    "outputs": {
                        "outfile": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/num_seqs.txt"  # noqa
                    },
                    "returnCode": 0,
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
                    "start": "2020-06-10T03:38:19.088Z",
                    "stderr": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout",  # noqa
                }
            ]
        },
        "end": "2020-06-10T03:43:56.709Z",
        "id": "ee30d68f-39d4-4fde-85c2-afdecce2bad3",
        "inputs": {
            "fq_count.fastq_file": "/global/cscratch1/sd/jaws/test/uploads/ekirton/NERSC/global/cfs/projectdirs/jaws/test/tiny.fastq"  # noqa
        },
        "labels": {
            "cromwell-workflow-id": "cromwell-ee30d68f-39d4-4fde-85c2-afdecce2bad3"
        },
        "outputs": {
            "fq_count.outfile": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/num_seqs.txt"  # noqa
        },
        "start": "2020-06-10T03:38:16.987Z",
        "status": "Succeeded",
        "submission": "2020-06-10T03:38:16.572Z",
        "submittedFiles": {
            "inputs": '{"fq_count.fastq_file":"/global/cscratch1/sd/jaws/test/uploads/ekirton/NERSC/global/cfs/projectdirs/jaws/test/tiny.fastq"}',  # noqa
            "labels": "{}",
            "options": "{\n\n}",
            "root": "",
            "workflow": 'workflow fq_count {\n    File fastq_file\n    call count_seqs { input: infile = fastq_file }\n    output {\n        File outfile = count_seqs.outfile\n    }\n}\n\ntask count_seqs {\n    File infile\n    command <<<\n        wc -l ${infile} | perl -ne \'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, " sequences\\n"} else {print STDERR "Invalid Fastq file\\n"}\' > num_seqs.txt\n    >>>\n    output {\n        File outfile = "num_seqs.txt"\n    }\n    runtime {\n        poolname: "test_small"\n        node: 1\n        nwpn: 1\n        mem: "10G"\n        time: "00:10:00"\n        shared: 0\n    }\n}\n',  # noqa
            "workflowUrl": "",
        },
        "workflowName": "fq_count",
        "workflowProcessingEvents": [
            {
                "cromwellId": "cromid-15c5397",
                "cromwellVersion": "47",
                "description": "PickedUp",
                "timestamp": "2020-06-10T03:38:16.985Z",
            },
            {
                "cromwellId": "cromid-15c5397",
                "cromwellVersion": "47",
                "description": "Finished",
                "timestamp": "2020-06-10T03:43:56.709Z",
            },
        ],
        "workflowRoot": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3",  # noqa
    },
    "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a": {  # METADATA FOR WORKFLOW_ID_EX2_MAIN
        "actualWorkflowLanguage": "WDL",
        "actualWorkflowLanguageVersion": "draft-2",
        "calls": {
            "main_workflow.goodbye": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-goodbye",  # noqa
                    "commandLine": 'echo "Goodbye World!"',
                    "end": "2020-09-10T23:07:47.221Z",
                    "executionEvents": [
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-10T23:07:47.221Z",
                            "startTime": "2020-09-10T23:07:46.572Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-10T23:07:36.343Z",
                            "startTime": "2020-09-10T23:07:36.332Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-10T23:07:36.331Z",
                            "startTime": "2020-09-10T23:07:35.804Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-10T23:07:35.804Z",
                            "startTime": "2020-09-10T23:07:35.803Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:36.332Z",
                            "startTime": "2020-09-10T23:07:36.331Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-10T23:07:46.572Z",
                            "startTime": "2020-09-10T23:07:36.343Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "World"},
                    "jobId": "5480",
                    "outputs": {"salutation": "Goodbye World!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "lr_jgicloud",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "condo_jgicloud",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:35.803Z",
                    "stderr": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-goodbye/execution/stderr",  # noqa
                    "stdout": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-goodbye/execution/stdout",  # noqa
                }
            ],
            "main_workflow.hello": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello",  # noqa
                    "commandLine": 'echo "Hello World!"',
                    "end": "2020-09-10T23:07:52.226Z",
                    "executionEvents": [
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-10T23:07:37.330Z",
                            "startTime": "2020-09-10T23:07:35.805Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:37.330Z",
                            "startTime": "2020-09-10T23:07:37.330Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-10T23:07:52.224Z",
                            "startTime": "2020-09-10T23:07:51.893Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-10T23:07:35.805Z",
                            "startTime": "2020-09-10T23:07:35.804Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-10T23:07:37.337Z",
                            "startTime": "2020-09-10T23:07:37.330Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-10T23:07:51.893Z",
                            "startTime": "2020-09-10T23:07:37.337Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "World"},
                    "jobId": "5481",
                    "outputs": {"salutation": "Hello World!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "lr_jgicloud",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "condo_jgicloud",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:35.804Z",
                    "stderr": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello/execution/stderr",  # noqa
                    "stdout": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello/execution/stdout",  # noqa
                }
            ],
            "main_workflow.hello_and_goodbye_1": [
                {
                    "attempt": 1,
                    "end": "2020-09-10T23:08:12.521Z",
                    "executionEvents": [
                        {
                            "description": "SubWorkflowPreparingState",
                            "endTime": "2020-09-10T23:07:35.807Z",
                            "startTime": "2020-09-10T23:07:35.805Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:35.805Z",
                            "startTime": "2020-09-10T23:07:35.804Z",
                        },
                        {
                            "description": "SubWorkflowRunningState",
                            "endTime": "2020-09-10T23:08:12.521Z",
                            "startTime": "2020-09-10T23:07:35.807Z",
                        },
                        {
                            "description": "SubWorkflowPendingState",
                            "endTime": "2020-09-10T23:07:35.804Z",
                            "startTime": "2020-09-10T23:07:35.804Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"hello_and_goodbye_input": "cruel world"},
                    "outputs": {
                        "goodbye_output": "Goodbye cruel world!",
                        "hello_output": "Hello cruel world!",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:35.804Z",
                    "subWorkflowId": "7408a4f1-bc85-49ba-8d5f-c886261ab6a0",
                }
            ],
            "main_workflow.hello_and_goodbye_2": [
                {
                    "attempt": 1,
                    "end": "2020-09-10T23:08:17.620Z",
                    "executionEvents": [
                        {
                            "description": "SubWorkflowPreparingState",
                            "endTime": "2020-09-10T23:07:35.808Z",
                            "startTime": "2020-09-10T23:07:35.805Z",
                        },
                        {
                            "description": "SubWorkflowPendingState",
                            "endTime": "2020-09-10T23:07:35.804Z",
                            "startTime": "2020-09-10T23:07:35.804Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:35.805Z",
                            "startTime": "2020-09-10T23:07:35.804Z",
                        },
                        {
                            "description": "SubWorkflowRunningState",
                            "endTime": "2020-09-10T23:08:17.620Z",
                            "startTime": "2020-09-10T23:07:35.808Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"hello_and_goodbye_input": "beautiful world"},
                    "outputs": {
                        "goodbye_output": "Goodbye beautiful world!",
                        "hello_output": "Hello beautiful world!",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:35.804Z",
                    "subWorkflowId": "89d86efc-dd04-48aa-a65f-21fb9d0c8be3",
                }
            ],
        },
        "end": "2020-09-10T23:08:19.665Z",
        "id": "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",
        "inputs": {},
        "labels": {
            "cromwell-workflow-id": "cromwell-74a0bf98-5bf3-4416-84bc-2fca6f4ed21a"
        },
        "metadataSource": "Unarchived",
        "outputs": {
            "main_workflow.main_output_1": "Hello cruel world!",
            "main_workflow.main_output_2": "Hello beautiful world!",
        },
        "start": "2020-09-10T23:07:33.678Z",
        "status": "Succeeded",
        "submission": "2020-09-10T23:07:33.637Z",
        "submittedFiles": {
            "imports": {
                "sub_wdl.wdl": 'task hello {\n  String addressee\n  command {\n    echo "Hello ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\ntask goodbye {\n  String addressee\n  command {\n    echo "Goodbye ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\nworkflow hello_and_goodbye {\n  String hello_and_goodbye_input\n\n  call hello {input: addressee = hello_and_goodbye_input }\n  call goodbye {input: addressee = hello_and_goodbye_input }\n\n  output {\n    String hello_output = hello.salutation\n    String goodbye_output = goodbye.salutation\n  }\n}\n'  # noqa
            },
            "inputs": "{}",
            "labels": "{}",
            "options": "{\n\n}",
            "root": "",
            "workflow": 'import "sub_wdl.wdl" as sub\n\nworkflow main_workflow {\n\n    call hello {input: addressee = "World" }\n    call goodbye {input: addressee = "World" }\n\n    call sub.hello_and_goodbye as hello_and_goodbye_1 { input: hello_and_goodbye_input = "cruel world" }\n    call sub.hello_and_goodbye as hello_and_goodbye_2 { input: hello_and_goodbye_input = "beautiful world" }\n\n    # call myTask { input: hello_and_goodbye.hello_output }\n\n    output {\n        String main_output_1 = hello_and_goodbye_1.hello_output\n        String main_output_2 = hello_and_goodbye_2.hello_output\n    }\n}\n\ntask hello {\n  String addressee\n  command {\n    echo "Hello ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\ntask goodbye {\n  String addressee\n  command {\n    echo "Goodbye ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n',  # noqa
            "workflowUrl": "",
        },
        "workflowName": "main_workflow",
        "workflowProcessingEvents": [
            {
                "cromwellId": "cromid-1366c4b",
                "cromwellVersion": "52",
                "description": "Finished",
                "timestamp": "2020-09-10T23:08:19.666Z",
            },
            {
                "cromwellId": "cromid-1366c4b",
                "cromwellVersion": "52",
                "description": "PickedUp",
                "timestamp": "2020-09-10T23:07:33.677Z",
            },
        ],
        "workflowRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",  # noqa
    },
    "89d86efc-dd04-48aa-a65f-21fb9d0c8be3": {  # METADATA FOR WORKFLOW_ID_EX2_SUB1
        "calls": {
            "hello_and_goodbye.goodbye": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_2/sub.hello_and_goodbye/89d86efc-dd04-48aa-a65f-21fb9d0c8be3/call-goodbye",  # noqa
                    "commandLine": 'echo "Goodbye beautiful world!"',
                    "end": "2020-09-10T23:08:12.221Z",
                    "executionEvents": [
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-10T23:08:11.345Z",
                            "startTime": "2020-09-10T23:07:39.339Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-10T23:07:39.339Z",
                            "startTime": "2020-09-10T23:07:39.330Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-10T23:08:12.221Z",
                            "startTime": "2020-09-10T23:08:11.345Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-10T23:07:39.330Z",
                            "startTime": "2020-09-10T23:07:37.841Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:39.330Z",
                            "startTime": "2020-09-10T23:07:39.330Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-10T23:07:37.841Z",
                            "startTime": "2020-09-10T23:07:37.840Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "beautiful world"},
                    "jobId": "5483",
                    "outputs": {"salutation": "Goodbye beautiful world!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "lr_jgicloud",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "condo_jgicloud",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:37.840Z",
                    "stderr": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_2/sub.hello_and_goodbye/89d86efc-dd04-48aa-a65f-21fb9d0c8be3/call-goodbye/execution/stderr",  # noqa
                    "stdout": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_2/sub.hello_and_goodbye/89d86efc-dd04-48aa-a65f-21fb9d0c8be3/call-goodbye/execution/stdout",  # noqa
                }
            ],
            "hello_and_goodbye.hello": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_2/sub.hello_and_goodbye/89d86efc-dd04-48aa-a65f-21fb9d0c8be3/call-hello",  # noqa
                    "commandLine": 'echo "Hello beautiful world!"',
                    "end": "2020-09-10T23:08:16.221Z",
                    "executionEvents": [
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-10T23:07:41.330Z",
                            "startTime": "2020-09-10T23:07:37.842Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-10T23:08:15.923Z",
                            "startTime": "2020-09-10T23:07:41.338Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:41.330Z",
                            "startTime": "2020-09-10T23:07:41.330Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-10T23:07:41.338Z",
                            "startTime": "2020-09-10T23:07:41.330Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-10T23:08:16.221Z",
                            "startTime": "2020-09-10T23:08:15.923Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-10T23:07:37.842Z",
                            "startTime": "2020-09-10T23:07:37.841Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "beautiful world"},
                    "jobId": "5485",
                    "outputs": {"salutation": "Hello beautiful world!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "lr_jgicloud",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "condo_jgicloud",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:37.841Z",
                    "stderr": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_2/sub.hello_and_goodbye/89d86efc-dd04-48aa-a65f-21fb9d0c8be3/call-hello/execution/stderr",  # noqa
                    "stdout": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_2/sub.hello_and_goodbye/89d86efc-dd04-48aa-a65f-21fb9d0c8be3/call-hello/execution/stdout",  # noqa
                }
            ],
        },
        "end": "2020-09-10T23:08:17.620Z",
        "id": "89d86efc-dd04-48aa-a65f-21fb9d0c8be3",
        "inputs": {"hello_and_goodbye_input": "beautiful world"},
        "metadataSource": "Unarchived",
        "outputs": {
            "hello_and_goodbye.goodbye_output": "Goodbye beautiful world!",
            "hello_and_goodbye.hello_output": "Hello beautiful world!",
        },
        "parentWorkflowId": "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",
        "rootWorkflowId": "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",
        "start": "2020-09-10T23:07:35.805Z",
        "status": "Succeeded",
        "workflowName": "hello_and_goodbye_2",
        "workflowRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",  # noqa
    },
    "7408a4f1-bc85-49ba-8d5f-c886261ab6a0": {  # METADATA FOR WORKFLOW_ID_EX2_SUB2
        "calls": {
            "hello_and_goodbye.goodbye": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_1/sub.hello_and_goodbye/7408a4f1-bc85-49ba-8d5f-c886261ab6a0/call-goodbye",  # noqa
                    "commandLine": 'echo "Goodbye cruel world!"',
                    "end": "2020-09-10T23:08:01.219Z",
                    "executionEvents": [
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-10T23:08:01.219Z",
                            "startTime": "2020-09-10T23:08:01.031Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-10T23:07:38.329Z",
                            "startTime": "2020-09-10T23:07:37.841Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-10T23:08:01.031Z",
                            "startTime": "2020-09-10T23:07:38.336Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-10T23:07:37.841Z",
                            "startTime": "2020-09-10T23:07:37.840Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:38.330Z",
                            "startTime": "2020-09-10T23:07:38.329Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-10T23:07:38.336Z",
                            "startTime": "2020-09-10T23:07:38.330Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "cruel world"},
                    "jobId": "5482",
                    "outputs": {"salutation": "Goodbye cruel world!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "lr_jgicloud",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "condo_jgicloud",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:37.840Z",
                    "stderr": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_1/sub.hello_and_goodbye/7408a4f1-bc85-49ba-8d5f-c886261ab6a0/call-goodbye/execution/stderr",  # noqa
                    "stdout": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_1/sub.hello_and_goodbye/7408a4f1-bc85-49ba-8d5f-c886261ab6a0/call-goodbye/execution/stdout",  # noqa
                }
            ],
            "hello_and_goodbye.hello": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_1/sub.hello_and_goodbye/7408a4f1-bc85-49ba-8d5f-c886261ab6a0/call-hello",  # noqa
                    "commandLine": 'echo "Hello cruel world!"',
                    "end": "2020-09-10T23:08:11.222Z",
                    "executionEvents": [
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-10T23:08:11.040Z",
                            "startTime": "2020-09-10T23:07:40.337Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-10T23:07:37.842Z",
                            "startTime": "2020-09-10T23:07:37.841Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-10T23:08:11.222Z",
                            "startTime": "2020-09-10T23:08:11.040Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-10T23:07:40.330Z",
                            "startTime": "2020-09-10T23:07:40.329Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-10T23:07:40.337Z",
                            "startTime": "2020-09-10T23:07:40.330Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-10T23:07:40.329Z",
                            "startTime": "2020-09-10T23:07:37.842Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "cruel world"},
                    "jobId": "5484",
                    "outputs": {"salutation": "Hello cruel world!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "lr_jgicloud",
                        "cluster": "jgi",
                        "constraint": "",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "false",
                        "maxRetries": "0",
                        "mem": "10G",
                        "node": "1",
                        "nwpn": "1",
                        "poolname": "test_small",
                        "qos": "condo_jgicloud",
                        "shared": "0",
                        "time": "00:10:00",
                    },
                    "shardIndex": -1,
                    "start": "2020-09-10T23:07:37.841Z",
                    "stderr": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_1/sub.hello_and_goodbye/7408a4f1-bc85-49ba-8d5f-c886261ab6a0/call-hello/execution/stderr",  # noqa
                    "stdout": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a/call-hello_and_goodbye_1/sub.hello_and_goodbye/7408a4f1-bc85-49ba-8d5f-c886261ab6a0/call-hello/execution/stdout",  # noqa
                }
            ],
        },
        "end": "2020-09-10T23:08:12.520Z",
        "id": "7408a4f1-bc85-49ba-8d5f-c886261ab6a0",
        "inputs": {"hello_and_goodbye_input": "cruel world"},
        "metadataSource": "Unarchived",
        "outputs": {
            "hello_and_goodbye.goodbye_output": "Goodbye cruel world!",
            "hello_and_goodbye.hello_output": "Hello cruel world!",
        },
        "parentWorkflowId": "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",
        "rootWorkflowId": "74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",
        "start": "2020-09-10T23:07:35.805Z",
        "status": "Succeeded",
        "workflowName": "hello_and_goodbye_1",
        "workflowRoot": "/global/scratch/jaws/jaws-dev/cromwell-executions/main_workflow/74a0bf98-5bf3-4416-84bc-2fca6f4ed21a",  # noqa
    },
    "dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15": {  # EX3 IS A FAILED RUN
        "actualWorkflowLanguage": "WDL",
        "actualWorkflowLanguageVersion": "draft-2",
        "calls": {
            "fq_count.count_seqs": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "callCaching": {
                        "allowResultReuse": False,
                        "effectiveCallCachingMode": "CallCachingOff",
                    },
                    "callRoot": "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs",  # noqa
                    "commandLine": 'wc -l /global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/inputs/-113930193/tiny.fastq | perl -ne \'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, " sequences\\n"} else {print STDERR "Invalid Fastq file\\n"}\' > num_seqs.txt',  # noqa
                    "end": "2020-09-11T22:30:34.735Z",
                    "executionEvents": [
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-09-11T22:30:33.918Z",
                            "startTime": "2020-09-11T22:30:33.853Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-09-11T22:30:34.736Z",
                            "startTime": "2020-09-11T22:30:34.257Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-09-11T22:30:33.566Z",
                            "startTime": "2020-09-11T22:30:33.554Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-09-11T22:30:34.257Z",
                            "startTime": "2020-09-11T22:30:33.918Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-09-11T22:30:33.853Z",
                            "startTime": "2020-09-11T22:30:33.845Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-09-11T22:30:33.845Z",
                            "startTime": "2020-09-11T22:30:33.566Z",
                        },
                    ],
                    "executionStatus": "Failed",
                    "failures": [
                        {
                            "causedBy": [],
                            "message": "Unable to start job. Check the stderr file for possible errors: /global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/execution/stderr.submit",  # noqa
                        }
                    ],
                    "inputs": {
                        "infile": "/global/cscratch1/sd/jaws_jtm/jaws-dev/uploads/ekirton/CORI/global/cfs/projectdirs/jaws/test/tiny.fastq"  # noqa
                    },
                    "retryableFailure": False,
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
                    "start": "2020-09-11T22:30:33.535Z",
                    "stderr": "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/execution/stdout",  # noqa
                }
            ]
        },
        "end": "2020-09-11T22:30:36.619Z",
        "failures": [
            {
                "causedBy": [
                    {
                        "causedBy": [],
                        "message": "Unable to start job. Check the stderr file for possible errors: /global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/execution/stderr.submit",  # noqa
                    }
                ],
                "message": "Workflow failed",
            }
        ],
        "id": "dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15",
        "inputs": {
            "fq_count.fastq_file": "/global/cscratch1/sd/jaws_jtm/jaws-dev/uploads/ekirton/CORI/global/cfs/projectdirs/jaws/test/tiny.fastq"  # noqa
        },
        "labels": {
            "cromwell-workflow-id": "cromwell-dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15"
        },
        "metadataSource": "Unarchived",
        "outputs": {},
        "start": "2020-09-11T22:30:30.640Z",
        "status": "Failed",
        "submission": "2020-09-11T22:30:29.607Z",
        "submittedFiles": {
            "inputs": '{"fq_count.fastq_file":"/global/cscratch1/sd/jaws_jtm/jaws-dev/uploads/ekirton/CORI/global/cfs/projectdirs/jaws/test/tiny.fastq"}',  # noqa
            "labels": "{}",
            "options": "{\n\n}",
            "root": "",
            "workflow": 'workflow fq_count {\n    File fastq_file\n    call count_seqs { input: infile = fastq_file }\n    output {\n        File outfile = count_seqs.outfile\n    }\n}\n\ntask count_seqs {\n    File infile\n    command <<<\n        wc -l ${infile} | perl -ne \'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, " sequences\\n"} else {print STDERR "Invalid Fastq file\\n"}\' > num_seqs.txt\n    >>>\n    output {\n        File outfile = "num_seqs.txt"\n    }\n    runtime {\n        poolname: "test_small"\n        node: 1\n        nwpn: 1\n        mem: "10G"\n        time: "00:10:00"\n    }\n}\n',  # noqa
            "workflowUrl": "",
        },
        "workflowName": "fq_count",
        "workflowProcessingEvents": [
            {
                "cromwellId": "cromid-41fef60",
                "cromwellVersion": "52",
                "description": "PickedUp",
                "timestamp": "2020-09-11T22:30:30.605Z",
            },
            {
                "cromwellId": "cromid-41fef60",
                "cromwellVersion": "52",
                "description": "Finished",
                "timestamp": "2020-09-11T22:30:36.620Z",
            },
        ],
        "workflowRoot": "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15",  # noqa
    },
}


def test_metadata():
    c = cromwell.Cromwell("localhost:8000")
    m = c.get_metadata(
        "ee30d68f-39d4-4fde-85c2-afdecce2bad3",
        METADATA["ee30d68f-39d4-4fde-85c2-afdecce2bad3"],
    )
    expectedWorkflowRoot = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # noqa
    workflowRoot = m.get("workflowRoot")
    assert workflowRoot == expectedWorkflowRoot


def test_task():
    EXAMPLE_TASK_NAME = "fq_count.count_seqs"
    EXAMPLE_CALLS = METADATA[WORKFLOW_ID_EX1]["calls"][EXAMPLE_TASK_NAME]
    task = cromwell.Task("localhost:8000", EXAMPLE_TASK_NAME, EXAMPLE_CALLS, METADATA)
    assert task.name == EXAMPLE_TASK_NAME
    assert task.is_subworkflow() is False
    assert task.calls[0]["attempt"] == 1
    assert task.calls[0]["jobId"] == "30"
    assert task.get("callRoot", 1) == "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs"  # noqa


def test_task_subworkflow():
    EXAMPLE_TASK_NAME = (
        "main_workflow.hello_and_goodbye_1"  # this task is a subworkflow
    )
    EXAMPLE_CALLS = METADATA[WORKFLOW_ID_EX2_MAIN]["calls"][EXAMPLE_TASK_NAME]
    task = cromwell.Task("localhost:8000", EXAMPLE_TASK_NAME, EXAMPLE_CALLS, METADATA)
    SUBWORKFLOW_JOB_IDS = (5482, 5484)
    SUBWORKFLOW_RUN_ID = "7408a4f1-bc85-49ba-8d5f-c886261ab6a0"
    for call in task.calls:
        if "jobId" in call:
            assert call["jobId"] in SUBWORKFLOW_JOB_IDS
        else:
            assert call["subWorkflowId"] == SUBWORKFLOW_RUN_ID


def test_get_all_metadata():
    """Given the uuid of a workflow with a subworkflow, we expect a dict with { uuid => metadata }"""
    c = cromwell.Cromwell("localhost:8000")
    result = c.get_all_metadata(WORKFLOW_ID_EX2_MAIN, METADATA)
    assert (
        bool(
            DeepDiff(
                result[WORKFLOW_ID_EX2_MAIN],
                METADATA[WORKFLOW_ID_EX2_MAIN],
                ignore_order=True,
            )
        )
        is False
    )
    assert (
        bool(
            DeepDiff(
                result[WORKFLOW_ID_EX2_SUB1],
                METADATA[WORKFLOW_ID_EX2_SUB1],
                ignore_order=True,
            )
        )
        is False
    )
    assert (
        bool(
            DeepDiff(
                result[WORKFLOW_ID_EX2_SUB2],
                METADATA[WORKFLOW_ID_EX2_SUB2],
                ignore_order=True,
            )
        )
        is False
    )


def test_get_errors():
    """Given workflow UUID, extract errors."""
    crom = cromwell.Cromwell("localhost:8000")
    metadata = crom.get_metadata(WORKFLOW_ID_EX3, METADATA[WORKFLOW_ID_EX3])
    status = metadata.execution_status()
    errors = metadata.errors()
    assert status["fq_count.count_seqs"] == "Failed"
    assert errors["fq_count.count_seqs"] == "Unable to start job. Check the stderr file for possible errors: /global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/fq_count/dfb4bc05-760d-4b0f-8a42-cc2fa3c78b15/call-count_seqs/execution/stderr.submit\nruntime:\n{'account': 'fungalp', 'cluster': 'cori', 'constraint': 'haswell', 'continueOnReturnCode': '0', 'cpu': '1', 'failOnStderr': 'false', 'maxRetries': '0', 'mem': '10G', 'node': '1', 'nwpn': '1', 'poolname': 'test_small', 'qos': 'genepool_special', 'shared': '0', 'time': '00:10:00'}\n" # noqa

    # test for no errors
    metadata = crom.get_metadata(WORKFLOW_ID_EX1, METADATA[WORKFLOW_ID_EX1])
    status = metadata.execution_status()
    errors = metadata.errors()
    assert errors == {}


def test_metadata_tasks():
    metadata = cromwell.Metadata(
        "localhost:8000", WORKFLOW_ID_EX2_MAIN, METADATA[WORKFLOW_ID_EX2_MAIN], METADATA
    )
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


def test_task_summary():
    metadata = cromwell.Metadata(
        "localhost:8000", WORKFLOW_ID_EX2_MAIN, METADATA[WORKFLOW_ID_EX2_MAIN], METADATA
    )
    result = metadata.task_summary()
    expected = [
        [WORKFLOW_ID_EX2_MAIN, "main_workflow.goodbye", 1, "5480"],
        [WORKFLOW_ID_EX2_MAIN, "main_workflow.hello", 1, "5481"],
        [WORKFLOW_ID_EX2_SUB2, "hello_and_goodbye.goodbye", 1, "5483"],
        [WORKFLOW_ID_EX2_SUB2, "hello_and_goodbye.hello", 1, "5485"],
        [WORKFLOW_ID_EX2_SUB1, "hello_and_goodbye.goodbye", 1, "5482"],
        [WORKFLOW_ID_EX2_SUB1, "hello_and_goodbye.hello", 1, "5484"],
    ]
    assert bool(DeepDiff(result, expected, ignore_order=True)) is False


def test_task_stdout():
    metadata = cromwell.Metadata("localhost:8000", WORKFLOW_ID_EX1, METADATA[WORKFLOW_ID_EX1], METADATA)
    expected_stderr = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr"  # noqa
    expected_stdout = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout"  # noqa
    assert len(metadata.tasks) == 1
    task = metadata.tasks[0]
    assert task.stderr() == expected_stderr
    assert task.stdout() == expected_stdout


def test_get_failed_task_runtime_attributes():
    expected_runtime = {
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
    }

    crom = cromwell.Cromwell("localhost:8000")
    metadata = crom.get_metadata(WORKFLOW_ID_EX3, METADATA[WORKFLOW_ID_EX3])
    for task in metadata.tasks:
        failures = task.failures()
        if failures:
            assert task.name == "fq_count.count_seqs"
            runtime = task.runtime()
            assert bool(DeepDiff(runtime, expected_runtime, ignore_order=True)) is False
