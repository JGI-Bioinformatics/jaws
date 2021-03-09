#!/usr/bin/env python
"""
These functions are to test the "testcases" from the "score_card" integration tests.
google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

This library of tests uses "fixtures" from conftest.py which should be located in the same directory. 
There is no need to import conftest.py as it is done automatically.

These test will cover issues that pertain to user debugging. 

Author: Jeff Froula <jlfroula@lbl.gov>
Updated: 02/23/21
"""

import pytest
import json
import time
from subprocess import Popen, PIPE
import submission_utils as util


#########################
###     Functions     ###
#########################
def test_see_result_failed_status(env,submit_bad_task):
    """ 
    When a user submits a WDL to site {param:site} and one of the tasks fail (for instance, due to a typo in a user-supplied command):
    the job needs to have result = failed (reflected by cmds: status, log, task-log, task-status)

    {
    "cromwell_run_id": "d0d045ca-6e9a-475c-87a9-96f76162f409",
    "download_task_id": "d7a54264-7628-11eb-8cfc-cd623f92e1c0",
    "id": 17028,
    "input_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    "input_site_id": "CORI",
    "output_dir": "/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/score-card-tests/o",
    "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    "result": "succeeded",
    }
    """

    # test status
    run_id = str(submit_bad_task['run_id'])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    data = json.loads(o) 
    assert data["result"] == 'failed'

def test_see_result_failed_task_status(env,submit_bad_task):
    """
    jaws run task-status 17028
    #TASK_NAME	ATTEMPT	CROMWELL_JOB_ID	STATUS_FROM	STATUS_TO	TIMESTAMP	REASON	STATUS_DETAIL
    fq_count.count_seqs	1	77273	running	failed	2021-02-23 22:45:04	failed with input file or command not found	The job has failed
    """
    # test task-status
    run_id = str(submit_bad_task['run_id'])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-status %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    assert 'failed with input file or command not found' in o

def test_see_result_failed_task_log(env,submit_bad_task):
    """
    jaws run task-log 17028
    #TASK_NAME	ATTEMPT	CROMWELL_JOB_ID	STATUS_FROM	STATUS_TO	TIMESTAMP	REASON
    fq_count.count_seqs	1	77273	created	ready	2021-02-23 22:44:10
    fq_count.count_seqs	1	77273	ready	queued	2021-02-23 22:44:10
    fq_count.count_seqs	1	77273	queued	pending	2021-02-23 22:44:12
    fq_count.count_seqs	1	77273	pending	running	2021-02-23 22:45:03
    fq_count.count_seqs	1	77273	running	failed	2021-02-23 22:45:04	failed with input file or command not found
    """

    # test task-log
    run_id = str(submit_bad_task['run_id'])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-log %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    assert 'failed with input file or command not found' in o

def test_see_result_failed_log(env,submit_bad_task):
    """
    jaws run log 17028
    #STATUS_FROM	STATUS_TO	TIMESTAMP	REASON
    created	uploading	2021-02-23 22:43:45	upload_task_id=967e420e-7628-11eb-8fff-01b9e52ec1df
    uploading	upload complete	2021-02-23 22:43:54
    upload complete	submitted	2021-02-23 22:44:06	cromwell_run_id=d0d045ca-6e9a-475c-87a9-96f76162f409
    submitted	queued	2021-02-23 22:44:18
    queued	running	2021-02-23 22:45:06
    running	succeeded	2021-02-23 22:45:22
    succeeded	downloading	2021-02-23 22:45:36	download_task_id=d7a54264-7628-11eb-8cfc-cd623f92e1c0
    downloading	download complete	2021-02-23 22:46:15
    """

    # test log
    run_id = str(submit_bad_task['run_id'])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run log %s | tail -n+2" % (env,run_id)
    (r,o,e) = util.run(cmd)

    a=[]
    for line in o.split("\n"):
        if line: 
            a.append(line.split()[1])
    
    assert 'failed' in a

def test_invalid_site(env,site):
    """
    jaws run submit WDLs/fq_count.wdl test-inputs/fq_count.json o smo

    SMO is not a valid Site ID.
    Available Sites:
    - {'max_ram_gb': '2048', 'site_id': 'CORI'}
    - {'max_ram_gb': '250', 'site_id': 'JGI'}

    """

    wdl="WDLs/fq_count.wdl"
    input_json="test-inputs/fq_count.json"
    outdir="./should-fail"

    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws run submit %s %s %s %s" % (env, wdl, input_json, outdir, "bogus")
    output = Popen(cmd, stdout=PIPE,
             stderr=PIPE, shell=True,
             universal_newlines=True)

    stdout,stderr=output.communicate()
    assert 'BOGUS is not a valid Site ID' in stdout

def test_invalid_docker_a(env,submit_bad_docker):
    """
    TESTCASE-33
    When user submits a wdl with a reference to a docker container that does not exist in the docker hub then:
    a) Job status should go to transition to failed
    b) error message should be available to user in the run's metadata
    """
    run_id = str(submit_bad_docker['run_id'])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-log %s" % (env,run_id)
    (r,o,e) = util.run(cmd)

    assert 'failed to run user command' in o

def test_invalid_docker_b(env,submit_bad_docker):
    """
    TESTCASE-33
    When user submits a wdl with a reference to a docker container that does not exist in the docker hub then:
    b) error message should be available to user in the run's metadata
    """
    run_id          = str(submit_bad_docker['run_id'])
    cromwell_run_id = submit_bad_docker['cromwell_run_id']
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run metadata %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    data = json.loads(o) 
    error_msg = data[cromwell_run_id]['failures'][0]['causedBy'][0]['message']
    assert 'docker not found' in error_msg

