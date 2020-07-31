from jaws_site import cromwell

CROMWELL_IDS = [
    "ee30d68f-39d4-4fde-85c2-afdecce2bad3",  # simple workflow without subworkflow
    "8595abd1-30c4-40e1-8de8-169e3df0209e",  # main workflow
    "ac7e42c1-456c-419f-a1e6-9e764bfc4af3",  # subworkflow of above
]

METADATA = {
    "ee30d68f-39d4-4fde-85c2-afdecce2bad3": {
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
                        "effectiveCallCachingMode": "CallCachingOff"
                    },
                    "callRoot": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs",  # noqa
                    "commandLine": "wc -l /global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/inputs/-230500445/tiny.fastq | perl -ne 'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, \" sequences\\n\"} else {print STDERR \"Invalid Fastq file\\n\"}' > num_seqs.txt",  # noqa
                    "end": "2020-06-10T03:43:54.794Z",
                    "executionEvents": [
                        {
                            "description": "RunningJob",
                            "endTime": "2020-06-10T03:43:54.673Z",
                            "startTime": "2020-06-10T03:38:19.896Z"
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-06-10T03:38:19.896Z",
                            "startTime": "2020-06-10T03:38:19.887Z"
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-06-10T03:38:19.887Z",
                            "startTime": "2020-06-10T03:38:19.887Z"
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-06-10T03:38:19.887Z",
                            "startTime": "2020-06-10T03:38:19.089Z"
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-06-10T03:43:54.794Z",
                            "startTime": "2020-06-10T03:43:54.673Z"
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-06-10T03:38:19.089Z",
                            "startTime": "2020-06-10T03:38:19.089Z"
                        }
                    ],
                    "executionStatus": "Done",
                    "inputs": {
                        "infile": "/global/cscratch1/sd/jaws/test/staging/ekirton/NERSC/global/cfs/projectdirs/jaws/test/tiny.fastq"  # noqa
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
                        "time": "00:10:00"
                    },
                    "shardIndex": -1,
                    "start": "2020-06-10T03:38:19.088Z",
                    "stderr": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3/call-count_seqs/execution/stdout"  # noqa
                }
            ]
        },
        "end": "2020-06-10T03:43:56.709Z",
        "id": "ee30d68f-39d4-4fde-85c2-afdecce2bad3",
        "inputs": {
            "fq_count.fastq_file": "/global/cscratch1/sd/jaws/test/staging/ekirton/NERSC/global/cfs/projectdirs/jaws/test/tiny.fastq"  # noqa
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
            "inputs": "{\"fq_count.fastq_file\":\"/global/cscratch1/sd/jaws/test/staging/ekirton/NERSC/global/cfs/projectdirs/jaws/test/tiny.fastq\"}",  # noqa
            "labels": "{}",
            "options": "{\n\n}",
            "root": "",
            "workflow": "workflow fq_count {\n    File fastq_file\n    call count_seqs { input: infile = fastq_file }\n    output {\n        File outfile = count_seqs.outfile\n    }\n}\n\ntask count_seqs {\n    File infile\n    command <<<\n        wc -l ${infile} | perl -ne 'if (/^\\s*(\\d+)/ and !($1%4)) {print $1/4, \" sequences\\n\"} else {print STDERR \"Invalid Fastq file\\n\"}' > num_seqs.txt\n    >>>\n    output {\n        File outfile = \"num_seqs.txt\"\n    }\n    runtime {\n        poolname: \"test_small\"\n        node: 1\n        nwpn: 1\n        mem: \"10G\"\n        time: \"00:10:00\"\n        shared: 0\n    }\n}\n",  # noqa
            "workflowUrl": ""
        },
        "workflowName": "fq_count",
        "workflowProcessingEvents": [
            {
                "cromwellId": "cromid-15c5397",
                "cromwellVersion": "47",
                "description": "PickedUp",
                "timestamp": "2020-06-10T03:38:16.985Z"
            },
            {
                "cromwellId": "cromid-15c5397",
                "cromwellVersion": "47",
                "description": "Finished",
                "timestamp": "2020-06-10T03:43:56.709Z"
            }
        ],
        "workflowRoot": "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # noqa
    },
    "8595abd1-30c4-40e1-8de8-169e3df0209e": {
        "actualWorkflowLanguage": "WDL",
        "actualWorkflowLanguageVersion": "draft-2",
        "calls": {
            "main_workflow.goodbye": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": True,
                        "effectiveCallCachingMode": "ReadAndWriteCache",
                        "hashes": {
                            "backend name": "6D3086C75F2DB761A86B2F982F10D384",
                            "command template": "ABFF2C1CC753BB668C97E9560BBDAF14",
                            "input": {
                                "String addressee": "362DC0726AA6A2A62E9F7773FD901DF1"
                            },
                            "input count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output expression": {
                                "String salutation": "0183144CF6617D5341681C6B2F756046"
                            },
                            "runtime attribute": {
                                "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
                                "docker": "N/A",
                                "failOnStderr": "68934A3E9455FA72420237EB05902327",
                            },
                        },
                        "hit": False,
                        "result": "Cache Miss",
                    },
                    "callRoot": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-goodbye",  # noqa
                    "commandLine": 'echo "Goodbye World!"',
                    "end": "2020-05-13T07:41:23.854Z",
                    "executionEvents": [
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-05-13T07:38:39.676Z",
                            "startTime": "2020-05-13T07:38:39.667Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-05-13T07:38:39.666Z",
                            "startTime": "2020-05-13T07:38:38.717Z",
                        },
                        {
                            "description": "UpdatingCallCache",
                            "endTime": "2020-05-13T07:41:22.974Z",
                            "startTime": "2020-05-13T07:41:22.461Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-05-13T07:38:39.667Z",
                            "startTime": "2020-05-13T07:38:39.666Z",
                        },
                        {
                            "description": "CallCacheReading",
                            "endTime": "2020-05-13T07:38:39.736Z",
                            "startTime": "2020-05-13T07:38:39.676Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-05-13T07:41:22.461Z",
                            "startTime": "2020-05-13T07:38:39.736Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-05-13T07:41:23.854Z",
                            "startTime": "2020-05-13T07:41:22.974Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-05-13T07:38:38.717Z",
                            "startTime": "2020-05-13T07:38:38.716Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "World"},
                    "jobId": "2597",
                    "outputs": {"salutation": "Goodbye World!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "fungalp",
                        "cluster": "cori",
                        "constraint": "haswell",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "False",
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
                    "start": "2020-05-13T07:38:38.716Z",
                    "stderr": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-goodbye/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-goodbye/execution/stdout",  # noqa
                }
            ],
            "main_workflow.hello": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": True,
                        "effectiveCallCachingMode": "ReadAndWriteCache",
                        "hashes": {
                            "backend name": "6D3086C75F2DB761A86B2F982F10D384",
                            "command template": "4EAADE3CD5D558C5A6CFA4FD101A1486",
                            "input": {
                                "String addressee": "362DC0726AA6A2A62E9F7773FD901DF1"
                            },
                            "input count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output expression": {
                                "String salutation": "0183144CF6617D5341681C6B2F756046"
                            },
                            "runtime attribute": {
                                "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
                                "docker": "N/A",
                                "failOnStderr": "68934A3E9455FA72420237EB05902327",
                            },
                        },
                        "hit": False,
                        "result": "Cache Miss",
                    },
                    "callRoot": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello",  # noqa
                    "commandLine": 'echo "Hello World!"',
                    "end": "2020-05-13T07:41:08.853Z",
                    "executionEvents": [
                        {
                            "description": "Pending",
                            "endTime": "2020-05-13T07:38:38.716Z",
                            "startTime": "2020-05-13T07:38:38.716Z",
                        },
                        {
                            "description": "CallCacheReading",
                            "endTime": "2020-05-13T07:38:39.707Z",
                            "startTime": "2020-05-13T07:38:39.702Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-05-13T07:41:08.853Z",
                            "startTime": "2020-05-13T07:41:07.967Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-05-13T07:41:06.281Z",
                            "startTime": "2020-05-13T07:38:39.707Z",
                        },
                        {
                            "description": "UpdatingCallCache",
                            "endTime": "2020-05-13T07:41:07.967Z",
                            "startTime": "2020-05-13T07:41:06.281Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-05-13T07:38:39.666Z",
                            "startTime": "2020-05-13T07:38:38.716Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-05-13T07:38:39.667Z",
                            "startTime": "2020-05-13T07:38:39.666Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-05-13T07:38:39.702Z",
                            "startTime": "2020-05-13T07:38:39.667Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "World"},
                    "jobId": "2596",
                    "outputs": {"salutation": "Hello World!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "fungalp",
                        "cluster": "cori",
                        "constraint": "haswell",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "False",
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
                    "start": "2020-05-13T07:38:38.716Z",
                    "stderr": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello/execution/stdout",  # noqa
                }
            ],
            "main_workflow.hello_and_goodbye": [
                {
                    "attempt": 1,
                    "end": "2020-05-13T07:53:52.742Z",
                    "executionEvents": [
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-05-13T07:38:38.739Z",
                            "startTime": "2020-05-13T07:38:38.735Z",
                        },
                        {
                            "description": "SubWorkflowPendingState",
                            "endTime": "2020-05-13T07:38:38.735Z",
                            "startTime": "2020-05-13T07:38:38.718Z",
                        },
                        {
                            "description": "SubWorkflowRunningState",
                            "endTime": "2020-05-13T07:53:52.738Z",
                            "startTime": "2020-05-13T07:38:38.755Z",
                        },
                        {
                            "description": "SubWorkflowPreparingState",
                            "endTime": "2020-05-13T07:38:38.755Z",
                            "startTime": "2020-05-13T07:38:38.739Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"hello_and_goodbye_input": "sub world"},
                    "outputs": {
                        "goodbye_output": "Goodbye sub world!",
                        "hello_output": "Hello sub world!",
                    },
                    "shardIndex": -1,
                    "start": "2020-05-13T07:38:38.715Z",
                    "subWorkflowId": "ac7e42c1-456c-419f-a1e6-9e764bfc4af3",
                }
            ],
        },
        "end": "2020-05-13T07:53:54.785Z",
        "id": "8595abd1-30c4-40e1-8de8-169e3df0209e",
        "inputs": {},
        "labels": {
            "cromwell-workflow-id": "cromwell-8595abd1-30c4-40e1-8de8-169e3df0209e"
        },
        "outputs": {"main_workflow.main_output": "Hello sub world!"},
        "start": "2020-05-13T07:38:36.463Z",
        "status": "Succeeded",
        "submission": "2020-05-13T07:38:36.363Z",
        "submittedFiles": {
            "imports": {
                "sub_wdl.wdl": 'task hello {\n  String addressee\n  command {\n    echo "Hello ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n      shared: 0\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\ntask goodbye {\n  String addressee\n  command {\n    echo "Goodbye ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n      shared: 0\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\nworkflow hello_and_goodbye {\n  String hello_and_goodbye_input\n\n  call hello {input: addressee = hello_and_goodbye_input }\n  call goodbye {input: addressee = hello_and_goodbye_input }\n\n  output {\n    String hello_output = hello.salutation\n    String goodbye_output = goodbye.salutation\n  }\n}\n'  # noqa
            },
            "inputs": "{}",
            "labels": "{}",
            "options": "{\n\n}",
            "root": "",
            "workflow": 'import "sub_wdl.wdl" as sub\n\nworkflow main_workflow {\n\n    call hello {input: addressee = "World" }\n    call goodbye {input: addressee = "World" }\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = "sub world" }\n\n    # call myTask { input: hello_and_goodbye.hello_output }\n\n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}\n\ntask hello {\n  String addressee\n  command {\n    echo "Hello ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n      shared: 0 \n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n\ntask goodbye {\n  String addressee\n  command {\n    echo "Goodbye ${addressee}!"\n  }\n  runtime {\n      poolname: "test_small"\n      node: 1\n      nwpn: 1\n      mem: "10G"\n      time: "00:10:00"\n      shared: 0\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n}\n',  # noqa
            "workflowUrl": "",
        },
        "workflowName": "main_workflow",
        "workflowProcessingEvents": [
            {
                "cromwellId": "cromid-275519a",
                "cromwellVersion": "47",
                "description": "Finished",
                "timestamp": "2020-05-13T07:53:54.785Z",
            },
            {
                "cromwellId": "cromid-275519a",
                "cromwellVersion": "47",
                "description": "PickedUp",
                "timestamp": "2020-05-13T07:38:36.462Z",
            },
        ],
        "workflowRoot": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e",  # noqa
    },
    "ac7e42c1-456c-419f-a1e6-9e764bfc4af3": {
        "calls": {
            "hello_and_goodbye.goodbye": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": True,
                        "effectiveCallCachingMode": "ReadAndWriteCache",
                        "hashes": {
                            "backend name": "6D3086C75F2DB761A86B2F982F10D384",
                            "command template": "ABFF2C1CC753BB668C97E9560BBDAF14",
                            "input": {
                                "String addressee": "B3EB934782A0B263A4A7759213CB3ADA"
                            },
                            "input count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output expression": {
                                "String salutation": "0183144CF6617D5341681C6B2F756046"
                            },
                            "runtime attribute": {
                                "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
                                "docker": "N/A",
                                "failOnStderr": "68934A3E9455FA72420237EB05902327",
                            },
                        },
                        "hit": False,
                        "result": "Cache Miss",
                    },
                    "callRoot": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello_and_goodbye/sub.hello_and_goodbye/ac7e42c1-456c-419f-a1e6-9e764bfc4af3/call-goodbye",  # noqa
                    "commandLine": 'echo "Goodbye sub world!"',
                    "end": "2020-05-13T07:53:29.856Z",
                    "executionEvents": [
                        {
                            "description": "Pending",
                            "endTime": "2020-05-13T07:38:40.789Z",
                            "startTime": "2020-05-13T07:38:40.786Z",
                        },
                        {
                            "description": "CallCacheReading",
                            "endTime": "2020-05-13T07:38:41.886Z",
                            "startTime": "2020-05-13T07:38:41.881Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-05-13T07:38:41.881Z",
                            "startTime": "2020-05-13T07:38:41.666Z",
                        },
                        {
                            "description": "UpdatingCallCache",
                            "endTime": "2020-05-13T07:53:28.969Z",
                            "startTime": "2020-05-13T07:53:27.787Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-05-13T07:38:41.666Z",
                            "startTime": "2020-05-13T07:38:41.666Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-05-13T07:53:29.856Z",
                            "startTime": "2020-05-13T07:53:28.969Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-05-13T07:38:41.666Z",
                            "startTime": "2020-05-13T07:38:40.789Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-05-13T07:53:27.787Z",
                            "startTime": "2020-05-13T07:38:41.886Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "sub world"},
                    "jobId": "2598",
                    "outputs": {"salutation": "Goodbye sub world!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "fungalp",
                        "cluster": "cori",
                        "constraint": "haswell",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "False",
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
                    "start": "2020-05-13T07:38:40.786Z",
                    "stderr": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello_and_goodbye/sub.hello_and_goodbye/ac7e42c1-456c-419f-a1e6-9e764bfc4af3/call-goodbye/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello_and_goodbye/sub.hello_and_goodbye/ac7e42c1-456c-419f-a1e6-9e764bfc4af3/call-goodbye/execution/stdout",  # noqa
                }
            ],
            "hello_and_goodbye.hello": [
                {
                    "attempt": 1,
                    "backend": "JTM",
                    "backendStatus": "Done",
                    "callCaching": {
                        "allowResultReuse": True,
                        "effectiveCallCachingMode": "ReadAndWriteCache",
                        "hashes": {
                            "backend name": "6D3086C75F2DB761A86B2F982F10D384",
                            "command template": "4EAADE3CD5D558C5A6CFA4FD101A1486",
                            "input": {
                                "String addressee": "B3EB934782A0B263A4A7759213CB3ADA"
                            },
                            "input count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output count": "C4CA4238A0B923820DCC509A6F75849B",
                            "output expression": {
                                "String salutation": "0183144CF6617D5341681C6B2F756046"
                            },
                            "runtime attribute": {
                                "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
                                "docker": "N/A",
                                "failOnStderr": "68934A3E9455FA72420237EB05902327",
                            },
                        },
                        "hit": False,
                        "result": "Cache Miss",
                    },
                    "callRoot": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello_and_goodbye/sub.hello_and_goodbye/ac7e42c1-456c-419f-a1e6-9e764bfc4af3/call-hello",  # noqa
                    "commandLine": 'echo "Hello sub world!"',
                    "end": "2020-05-13T07:53:50.861Z",
                    "executionEvents": [
                        {
                            "description": "UpdatingCallCache",
                            "endTime": "2020-05-13T07:53:49.972Z",
                            "startTime": "2020-05-13T07:53:48.119Z",
                        },
                        {
                            "description": "CallCacheReading",
                            "endTime": "2020-05-13T07:38:41.891Z",
                            "startTime": "2020-05-13T07:38:41.886Z",
                        },
                        {
                            "description": "PreparingJob",
                            "endTime": "2020-05-13T07:38:41.886Z",
                            "startTime": "2020-05-13T07:38:41.666Z",
                        },
                        {
                            "description": "UpdatingJobStore",
                            "endTime": "2020-05-13T07:53:50.861Z",
                            "startTime": "2020-05-13T07:53:49.972Z",
                        },
                        {
                            "description": "RunningJob",
                            "endTime": "2020-05-13T07:53:48.119Z",
                            "startTime": "2020-05-13T07:38:41.891Z",
                        },
                        {
                            "description": "RequestingExecutionToken",
                            "endTime": "2020-05-13T07:38:41.666Z",
                            "startTime": "2020-05-13T07:38:40.786Z",
                        },
                        {
                            "description": "Pending",
                            "endTime": "2020-05-13T07:38:40.786Z",
                            "startTime": "2020-05-13T07:38:40.786Z",
                        },
                        {
                            "description": "WaitingForValueStore",
                            "endTime": "2020-05-13T07:38:41.666Z",
                            "startTime": "2020-05-13T07:38:41.666Z",
                        },
                    ],
                    "executionStatus": "Done",
                    "inputs": {"addressee": "sub world"},
                    "jobId": "2599",
                    "outputs": {"salutation": "Hello sub world!"},
                    "returnCode": 0,
                    "runtimeAttributes": {
                        "account": "fungalp",
                        "cluster": "cori",
                        "constraint": "haswell",
                        "continueOnReturnCode": "0",
                        "cpu": "1",
                        "failOnStderr": "False",
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
                    "start": "2020-05-13T07:38:40.786Z",
                    "stderr": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello_and_goodbye/sub.hello_and_goodbye/ac7e42c1-456c-419f-a1e6-9e764bfc4af3/call-hello/execution/stderr",  # noqa
                    "stdout": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e/call-hello_and_goodbye/sub.hello_and_goodbye/ac7e42c1-456c-419f-a1e6-9e764bfc4af3/call-hello/execution/stdout",  # noqa
                }
            ],
        },
        "end": "2020-05-13T07:53:52.738Z",
        "id": "ac7e42c1-456c-419f-a1e6-9e764bfc4af3",
        "inputs": {"hello_and_goodbye_input": "sub world"},
        "outputs": {
            "hello_and_goodbye.goodbye_output": "Goodbye sub world!",
            "hello_and_goodbye.hello_output": "Hello sub world!",
        },
        "parentWorkflowId": "8595abd1-30c4-40e1-8de8-169e3df0209e",
        "rootWorkflowId": "8595abd1-30c4-40e1-8de8-169e3df0209e",
        "start": "2020-05-13T07:38:38.738Z",
        "status": "Succeeded",
        "workflowName": "hello_and_goodbye",
        "workflowRoot": "/global/cscratch1/sd/jaws_jtm/dev/cromwell-executions/main_workflow/8595abd1-30c4-40e1-8de8-169e3df0209e",  # noqa
    },
}


def test_metadata():
    c = cromwell.Cromwell("localhost:8000")
    m = c.get_metadata(
        "ee30d68f-39d4-4fde-85c2-afdecce2bad3",
        METADATA["ee30d68f-39d4-4fde-85c2-afdecce2bad3"]
    )
    expectedWorkflowRoot = "/global/cscratch1/sd/jaws/test/cromwell-executions/fq_count/ee30d68f-39d4-4fde-85c2-afdecce2bad3"  # noqa
    workflowRoot = m.get("workflowRoot")
    assert workflowRoot == expectedWorkflowRoot
