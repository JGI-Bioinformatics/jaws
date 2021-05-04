#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.

import os
import pytest
import json
import time
import submission_utils as util

tmp_wdl = "pow23.wdl"
tmp_readme = "pow23.md"
released_wdl_catalog_name="fq_count"
released_wdl_catalog_wdl="fq_count.wdl"
changed_readme = "fq_count_new.md"
check_tries=100
check_sleep=30


#########################
###     Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#
def test_jaws_info(env):
    """ tests that there is a valid output for jaws info. Name should be dev,staging, or prod and version should have some value.
    {
    "docs_url": "https://jaws-docs.readthedocs.io/en/latest/",
        "name": "prod",
        "version": "2.1"
    } """

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws info" % (env)
    (r,o,e) = util.run(cmd)

    data = json.loads(o)

    # do we have an acceptable name
    assert data["name"] in ["prod","staging","dev"]
    assert data["version"] is not None

def test_jaws_status(env):
    """ tests that the jaws status is working. We don't care if some services are down.
        Just test that all below services are shown, regardless of status.
    {
        "CORI-Cromwell": "UP",
        "CORI-RMQ": "UP",
        "CORI-Site": "UP",
        "JAWS-Central": "UP",
        "JGI-Cromwell": "UP",
        "JGI-RMQ": "UP",
        "JGI-Site": "UP"
    }
    """

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws health" % (env)
    (r,o,e) = util.run(cmd)
    data = json.loads(o)

    actual_keys = list(data.keys())
    required_keys = ["CORI-Cromwell","CORI-RMQ","CORI-Site","JAWS-Central","JGI-Cromwell","JGI-RMQ","JGI-Site"]

    for k in required_keys:
        assert k in actual_keys


def test_jaws_run_queue(env, submit_fq_count_wdl):
    """ tests that the jaws queue command has the run id in the stdout.
    [
    {
        "id": 5759,
        "input_site_id": "CORI",
        "json_file": "/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/deployment-tests/test-inputs/fq_count.json",
        "result": null,
        "site_id": "CORI",
        "status": "submitted",
        "status_detail": "The run has been submitted to Cromwell and tasks should start to queue within moments",
        "submitted": "2021-05-04 01:44:36",
        "tag": null,
        "updated": "2021-05-04 01:45:03",
        "wdl_file": "/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/deployment-tests/WDLs/fq_count.wdl"
    }
    ]
    """

    run_id = str(submit_fq_count_wdl['run_id'])

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws queue | grep '\"id\":' | awk '{print $2}' | tr -d ','" % (env)
    (r,o,e) = util.run(cmd)
    ids=o.split()
    assert run_id in ids


def test_jaws_run_history(env, submit_fq_count_wdl):
    """ tests that the jaws history command has the run id in the stdout."""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws history | grep '\"id\": %s' | awk '{print $2}' | tr -d ','" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert o

def test_jaws_wdl_metadata(env, submit_fq_count_wdl):
    """Check that a jaws metadata returns workflowRoot has a value"""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws metadata %s | grep workflowRoot | awk '{print $2}' | tr -d '\"'" % (env,run_id)
    (r,o,e) = util.run(cmd)

    assert o 

def test_jaws_wdl_errors(env, submit_fq_count_wdl):
    """Check that a jaws metadata returns workflowRoot has a value"""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws errors %s" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # we don't have any errors in our wdl so the errors command should't return anything.
    # We don't have a good test for this command here, but we only test that the return code is 0.
    assert r == 0
        
def test_jaws_wdl_task_status(env, submit_fq_count_wdl):
    """Check that jaws task-status returns something like this:
     fq_count.count_seqs 1   25177   running success 2021-01-13 12:37:45     The job completed successfully

     It should have only one line.
    """
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)
    #time.sleep(120)  # wait an additional amount of time to make sure everything is updated

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws task-status %s" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # make sure there is only 1 output line (not including header)
    a=o.split("\n")
    a = list(filter(None,a)) # remove empty elements
    assert len(a) == 2

    # output line should have this string
    assert 'fq_count.count_seqs' in a[1]

def test_jaws_wdl_log(env, submit_fq_count_wdl):
    """Check that the first line of jaws log returns something like this:
       created uploading 2021-04-06 02:56:49 upload_task_id=bbbc09c2-9683-11eb-955a-752ba7b88ebe
    """
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)
    #time.sleep(120)  # wait an additional amount of time to make sure everything is updated

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws log %s " % (env,run_id)
    (r,o,e) = util.run(cmd)
    a=o.split("\n")
    a = list(filter(None,a)) # remove empty elements

    assert 'created' in a[1] and 'uploading' in a[1]

    # there should be 4 columns of output (but time is split into two so actually 5)
    assert len(a[1].split()) == 5

    # there should be 9 output rows (including header)
    assert len(a) == 9


def test_jaws_wdl_task_log(env, submit_fq_count_wdl):
    """Check that the first line of jaws task-log returns something like this:
       f0b7fd65-1620-4765-8b62-55d0bec74a8d  fq_count.count_seqs 1 16948  running failed 2021-04-06 03:46:26
    """
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)
    #time.sleep(120)  # wait an aditional amount of time to make sure everything is updated

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws task-log %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    a=o.split("\n")
    a = list(filter(None,a)) # remove empty elements

    assert 'fq_count.count_seqs' in a[1]

    # there should be 8 columns of output 
    assert len(a[1].split()) == 8

    # there should be 6 output rows (including header)
    assert len(a) == 6


def mtest_wfcopy(env,dir,submit_fq_count_wdl):
    """ 
    Check that wfcopy works and can flatten the directory

    i.e. should see something like MyCopy/count_seqs/num_seqs.txt
    """
    # remove any old copy of the temp dir if exists
    mycopy = "./Mycopy"
    if os.path.exists(mycopy):
        cmd = "rm -rf %s" % mycopy
        (r,o,e) = util.run(cmd)
        assert not r

    run_id = str(submit_fq_count_wdl['run_id'])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status --verbose %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    data = json.loads(o)
    outdir = data['output_dir']

    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wfcopy --flatten %s %s" % (env,outdir,mycopy)
    (r,o,e) = util.run(cmd)
    assert not r

    # does this path exist
    assert os.path.exists(os.path.join(mycopy,"count_seqs/num_seqs.txt"))


def test_jaws_queue_site_filter(env, site, submit_fq_count_wdl):
    """
    jaws queue --site [CORI, JGI, CASCADE]
    """
    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws queue --site %s" % (env, site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["site_id"].lower() == site
    else:
        pytest.exit(f"no runs were found in the queue for site: {site}")


def test_jaws_history_site_filter(env, site, submit_fq_count_wdl):
    """
    jaws history --site [CORI, JGI, CASCADE]
    """
    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --site %s" % (env, site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["site_id"].lower() == site
    else:
        pytest.exit(f"no runs were found in the history for site: {site}")


def test_jaws_history_result_filter_succeeded(env,submit_fq_count_wdl):
    """
    jaws history --result [succeeded, failed]
    Checking the output only with "succeeded" and "failed"
    """

    run_id = submit_fq_count_wdl['run_id']
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --result succeeded" % (env)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["result"] == 'succeeded'
    else:
        pytest.exit(f"no runs were found in the history that have result: succeeded")


def test_jaws_history_result_filter_failed(env,submit_bad_task):
    """
    jaws history --result failed
    Checking the output only with "succeeded" and "failed"
    """
    run_id = submit_bad_task['run_id']
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --result failed" % (env)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["result"] == 'failed'
    else:
        pytest.exit(f"no runs were found in the history that have result: failed")

