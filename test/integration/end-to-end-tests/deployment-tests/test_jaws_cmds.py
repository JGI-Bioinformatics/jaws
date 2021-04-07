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

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status" % (env)
    (r,o,e) = util.run(cmd)
    data = json.loads(o)

    actual_keys = list(data.keys())
    required_keys = ["CORI-Cromwell","CORI-RMQ","CORI-Site","JAWS-Central","JGI-Cromwell","JGI-RMQ","JGI-Site"]

    for k in required_keys:
        assert k in actual_keys


def test_jaws_run_queue(env, submit_fq_count_wdl):
    """ tests that the jaws run queue command has the run id in the stdout."""

    data = submit_fq_count_wdl
    run_id = str(data['run_id'])

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run queue | grep '\"id\":' | awk '{print $2}' | tr -d ','" % (env)
    (r,o,e) = util.run(cmd)
    ids=o.split()
    assert run_id in ids


def test_jaws_run_history(env, submit_fq_count_wdl):
    """ tests that the jaws run history command has the run id in the stdout."""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run history | grep '\"id\": %s' | awk '{print $2}' | tr -d ','" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert o

def test_jaws_wdl_update_readme(env,dir):
    """This tests for multiple things 
        1) adding a WDL to the catalog, 
        2) updating it's readme,  
        3) updating the WDL. 
        4) deleted the WDL. 
    """
    ori_wdl=os.path.join(dir,"WDLs",released_wdl_catalog_wdl)
    ori_readme = os.path.join(dir,"test-inputs/fq_count.md")
    changed_wdl=os.path.join(dir,"WDLs/fq_count_changed.wdl")

    with open(changed_readme,"w") as f:
        f.write('this readme has been changed')

    # make sure wdl doesn't already exist in the catalog
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list" % (env)
    (r,o,e) = util.run(cmd)
    if 'uniq_version' in o:
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl delete %s uniq_version" % (env,released_wdl_catalog_name)
        (r,o,e) = util.run(cmd)
        assert not r

    # add to catalog
    # this command returns a rc > 0 even upon success, so I can't test for rc, but use "result" from stderr.
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl add %s uniq_version %s %s" % (env,released_wdl_catalog_name,ori_wdl,ori_readme)
    (r,o,e) = util.run(cmd)
    data = json.loads(e)
    assert data['result'] == 'OK'

    # update-doc
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-doc %s uniq_version %s" % (env,released_wdl_catalog_name,changed_readme)
    (r,o,e) = util.run(cmd)
    data = json.loads(o)
    assert data['result'] == 'OK'
    
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl about %s uniq_version" % (env,released_wdl_catalog_name)
    (r,o,e) = util.run(cmd)
    assert 'changed' in o

    # update-wdl
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-wdl %s uniq_version %s" % (env,released_wdl_catalog_name,changed_wdl)
    (r,o,e) = util.run(cmd)
    data = json.loads(o)
    assert data['result'] == 'OK'
    
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl get %s uniq_version " % (env,released_wdl_catalog_name)
    (r,o,e) = util.run(cmd)
    assert 'fq_count_changed' in o

    # delete WDL
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl delete %s uniq_version" % (env,released_wdl_catalog_name)
    (r,o,e) = util.run(cmd)
    assert not r

def test_jaws_wdl_update_released_wdl(env,dir):
    """You should not be able to change a WDL that has been released"""
    changed_wdl=os.path.join(dir,"WDLs/fq_count_changed.wdl")

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-wdl %s 1.0.0 %s" % (env,released_wdl_catalog_name,changed_wdl)
    (r,o,e) = util.run(cmd)
    assert "Action not allowed" in e

def test_jaws_wdl_versions(env):
    """Test that we can get the version for a given WDL
    "fq_count:1.0.0": {
        "created": "2020-11-04T22:24:11Z",
        "last_updated": "2020-11-04T22:24:11Z",
        "name": "fq_count",
        "owner": "jfroula",
        "production_release": "no",
        "version": "1.0.0"
        }
    """
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl versions %s | grep version | awk '{print $2}' " % (env,released_wdl_catalog_name)
    (r,o,e) = util.run(cmd)
    assert o.strip() == '\"1.0.0\"'

def test_jaws_wdl_metadata(env, submit_fq_count_wdl):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run metadata %s | grep workflowRoot | awk '{print $2}' | tr -d '\"'" % (env,run_id)
    (r,o,e) = util.run(cmd)

    assert o 

def test_jaws_wdl_errors(env, submit_fq_count_wdl):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run errors %s" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # we don't have any errors in our wdl so the errors command should't return anything.
    # We don't have a good test for this command here, but we only test that the return code is 0.
    assert r == 0
        
def test_jaws_wdl_task_status(env, submit_fq_count_wdl):
    """Check that jaws run task-status returns something like this:
     fq_count.count_seqs 1   25177   running success 2021-01-13 12:37:45     The job completed successfully

     It should have only one line.
    """
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)
    #time.sleep(120)  # wait an additional amount of time to make sure everything is updated

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-status %s" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # make sure there is only 1 output line (not including header)
    a=o.split("\n")
    a = list(filter(None,a)) # remove empty elements
    assert len(a) == 2

    # output line should have this string
    assert 'fq_count.count_seqs' in a[1]

def test_jaws_wdl_log(env, submit_fq_count_wdl):
    """Check that the first line of jaws run log returns something like this:
       created uploading 2021-04-06 02:56:49 upload_task_id=bbbc09c2-9683-11eb-955a-752ba7b88ebe
    """
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)
    #time.sleep(120)  # wait an additional amount of time to make sure everything is updated

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run log %s " % (env,run_id)
    (r,o,e) = util.run(cmd)
    a=o.split("\n")
    a = list(filter(None,a)) # remove empty elements

    assert 'created' in a[1] and 'uploading' in a[1]

    # there should be 4 columns of output (but time is split into two so actually 5)
    assert len(a[1].split()) == 5

    # there should be 9 output rows (including header)
    assert len(a) == 9


def mtest_jaws_wdl_task_log(env, submit_fq_count_wdl):
    """Check that the first line of jaws run task-log returns something like this:
       f0b7fd65-1620-4765-8b62-55d0bec74a8d  fq_count.count_seqs 1 16948  running failed 2021-04-06 03:46:26
    """
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(run_id,env,check_tries,check_sleep)
    #time.sleep(120)  # wait an aditional amount of time to make sure everything is updated

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-log %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    a=o.split("\n")
    a = list(filter(None,a)) # remove empty elements

    assert 'fq_count.count_seqs' in a[1]

    # there should be 8 columns of output 
    assert len(a[1].split()) == 8

    # there should be 6 output rows (including header)
    assert len(a) == 6


def test_wfcopy(env,dir,submit_fq_count_wdl):
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
    outdir = submit_fq_count_wdl['output_dir']
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wfcopy --flatten %s %s" % (env,outdir,mycopy)
    (r,o,e) = util.run(cmd)
    assert not r

    # does this path exist
    assert os.path.exists(os.path.join(mycopy,"count_seqs/num_seqs.txt"))

