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

import json
import time
from subprocess import Popen, PIPE
import submission_utils as util


#####################
#     Functions     #
#####################
def test_should_fail_status(submit_bad_task):
    """
    When a user submits a WDL to site {param:site} and one of the tasks
    fail (for instance, due to a typo in a user-supplied command):
    the job needs to have result = failed (reflected by cmds: status, log, task-log, task-status)

    {
    "cromwell_id": "d0d045ca-6e9a-475c-87a9-96f76162f409",
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
    run_id = str(submit_bad_task["run_id"])
    cmd = "jaws status --verbose %s" % (run_id)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    assert data["result"] == "failed", "jaws-status should say run failed"


def test_should_fail_task_status(submit_bad_task):
    """
    jaws task-status 17028
    #TASK_NAME  CROMWELL_JOB_ID STATUS       TIMESTAMP       REASON
    fq_count.count_seqs 77273   failed  2021-02-23 22:45:04     failed with input file or command not found
    """ # noqa
    # test task-status
    id = str(submit_bad_task["run_id"])
    cmd = "jaws task-status %s" % (id)
    (r, o, e) = util.run(cmd)
    assert "command not found" in o.replace("\n", " ")


def test_should_fail_task_log(submit_bad_task):
    """
    jaws task-log 17028
    #TASK_NAME  CROMWELL_JOB_ID STATUS_FROM     STATUS_TO       TIMESTAMP       REASON
    fq_count.count_seqs 77273   created ready   2021-02-23 22:44:10
    fq_count.count_seqs 77273   ready   queued  2021-02-23 22:44:10
    fq_count.count_seqs 77273   queued  pending 2021-02-23 22:44:12
    fq_count.count_seqs 77273   pending running 2021-02-23 22:45:03
    fq_count.count_seqs 77273   running failed  2021-02-23 22:45:04     failed with input file or command not found
    """ # noqa

    # test task-log
    id = str(submit_bad_task["run_id"])
    cmd = "jaws task-log %s" % (id)
    (r, o, e) = util.run(cmd)
    assert "command not found" in o.replace("\n", " ")


def test_should_fail_log(submit_bad_task):
    """
    jaws log 17028
    #STATUS_FROM        STATUS_TO       TIMESTAMP       REASON
    created     uploading       2021-02-23 22:43:45     upload_task_id=967e420e-7628-11eb-8fff-01b9e52ec1df
    uploading   upload complete 2021-02-23 22:43:54
    upload complete     submitted       2021-02-23 22:44:06     cromwell_id=d0d045ca-6e9a-475c-87a9-96f76162f409
    submitted   queued  2021-02-23 22:44:18
    queued      running 2021-02-23 22:45:06
    running     succeeded       2021-02-23 22:45:22
    succeeded   downloading     2021-02-23 22:45:36     download_task_id=d7a54264-7628-11eb-8cfc-cd623f92e1c0
    downloading download complete       2021-02-23 22:46:15
    """

    # test log
    id = str(submit_bad_task["run_id"])
    cmd = "jaws log %s | tail -n+2" % (id)
    (r, o, e) = util.run(cmd)

    a = []
    for line in o.split("\n"):
        if line:
            a.append(line.split()[1])

    assert "failed" in a, "jaws-log should say run failed"


def test_invalid_site(site):
    """
    jaws submit --no-cache --quiet WDLs/fq_count.wdl test-inputs/fq_count.json o smo

    SMO is not a valid Site ID.
    Available Sites:
    - {'max_ram_gb': '2048', 'site_id': 'CORI'}
    - {'max_ram_gb': '250', 'site_id': 'JGI'}

    example error output
    {'detail': {'error': 'Unknown Site ID; "BOSUG" is not one of our sites'}, 'status': 404, 'title': 'Not Found', 'type': 'about:blank'}
    """ # noqa

    wdl = "WDLs/fq_count.wdl"
    input_json = "test-inputs/fq_count.json"

    cmd = "jaws submit --quiet --no-cache %s %s %s" % (
        wdl,
        input_json,
        "bogus",
    )
    output = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, universal_newlines=True)

    stdout, stderr = output.communicate()
    assert (
        'Unknown Site ID; "BOGUS"' in stderr
    ), "bogus is not an acceptable site, run should fail"


def test_invalid_docker_a(submit_bad_docker):
    """
    TESTCASE-33
    When user submits a wdl with a reference to a docker container that does not exist in the docker hub then:
    a) Job status should go to transition to failed
    """
    id = str(submit_bad_docker["run_id"])
    cmd = "jaws log %s" % (id)
    (r, o, e) = util.run(cmd)

    assert "failed" in o, "jaws-log should say run failed"


def test_invalid_docker_b(site, submit_bad_docker):
    """
    TESTCASE-33
    When user submits a wdl with a reference to a docker container that does not exist in the docker hub then:
    b) error message should be available to user in the run's metadata
    """
    # get cromwell id from status
    id = str(submit_bad_docker["run_id"])
    cmd = "jaws status --verbose %s" % (id)
    (r, o, e) = util.run(cmd)

    # check the metadata
    cmd = "jaws errors %s" % (id)
    (r, o, e) = util.run(cmd)
    if site.lower() == 'cori':
        assert (
            "Invalid container name or failed to pull container" in o
        ), "There should be a message saying docker was not found"
    elif site.lower() == 'jgi' or site.lower() == 'tahoma':
        assert (
            "Failed to pull" in o
        ), "There should be a message saying docker was not found"
    else:
        assert 0, f"Expected site to be cori or jgi but found {site}"


def test_bad_sub_workflow_error_msg(submit_bad_sub_task):
    """
    TESTCASE-42
    When user submits a wdl with a subworkflow that has a command error then
    the errors command should display error message
    """
    id = str(submit_bad_sub_task["run_id"])
    cmd = "jaws errors %s" % (id)
    (r, o, e) = util.run(cmd)

    assert "echoooo: command not found" in o, "sub workflow command error should appear in errors"


def test_timeout(dir, site):
    """
    TESTCASE-44
    When user submits a wdl and a timeout occurs, the
    timeout message should appear in the output from the errors command
    """
    WDL = "/WDLs/timeout.wdl"
    INP = "/test-inputs/timeout.json"
    check_sleep = 30
    check_tries = 50

    wdl = dir + WDL
    input_json = dir + INP

    run_id = util.submit_wdl(wdl, input_json, site)["run_id"]
    util.wait_for_run(run_id, check_tries, check_sleep)

    time.sleep(60)

    # get the errors from JAWS for that run
    cmd = "jaws errors %s" % (run_id)
    r, o, e = util.run(cmd)

    # do the check!
    fail_msg = "error. Keyword absent: \"timeout\" (%s)" % run_id
    assert "failed with timeout" in o, fail_msg


def test_bad_ref_dir(dir, site):
    """
    When user submits a wdl with a bad path to /refdata like:
    /refdata/i_dont_exist
    We should get a user friendly error message like:

      No such file or directory
    """
    WDL = "/WDLs/bad_ref.wdl"
    INP = "/test-inputs/bad_ref.json"
    check_sleep = 30
    check_tries = 50

    wdl = dir + WDL
    input_json = dir + INP

    run_id = util.submit_wdl(wdl, input_json, site)["run_id"]
    util.wait_for_run(run_id, check_tries, check_sleep)

    # get the errors from JAWS for that run
    cmd = "jaws errors %s" % (run_id)
    r, o, e = util.run(cmd)

    # do the check!
    fail_msg = "Error should say: \"No such file or directory\""
    assert "No such file or directory" in o, fail_msg
