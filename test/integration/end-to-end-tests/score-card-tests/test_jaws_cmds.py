#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.

import pytest
import json
import time
import parsing_functions as pf

tmp_wdl = "pow23.wdl"
tmp_readme = "pow23.md"
wdl_catalog_name="tmp_wdl_catalog_name"
released_wdl_catalog_name="fq_count"
check_tries=100
check_sleep=30


@pytest.fixture
def command_line_args(submission_info, env, site):
    submission_info["env"] = env
    submission_info["site"] = site


#########################
###     Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#
def wait_for_run(env,run_id):
    """ Wait for all the runs in run_ids list to finish."""
    tries = 1 
    while tries <= check_tries:
        # check whether the run has finished every 60 seconds
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env,run_id)
        (o,e,r) = pf.submit_cmd(cmd)
        if r > 0:
            pytest.exit("stderr: %s" % e)

        status_output = json.loads(o)
        run_status = status_output["status"]

        if run_status == "download complete":
            return

        tries += 1
        time.sleep(check_sleep)

def test_jaws_info(env):
    """ tests that there is a valid output for jaws info. Name should be dev,staging, or prod and version should have some value.
    {
    "docs_url": "https://jaws-docs.readthedocs.io/en/latest/",
        "name": "prod",
        "version": "2.1"
    } """

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws info" % (env)
    (o,e,r) = pf.submit_cmd(cmd)

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
    (o,e,r) = pf.submit_cmd(cmd)
    data = json.loads(o)

    actual_keys = list(data.keys())
    required_keys = ["CORI-Cromwell","CORI-RMQ","CORI-Site","JAWS-Central","JGI-Cromwell","JGI-RMQ","JGI-Site"]

    for k in required_keys:
        assert k in actual_keys


def test_jaws_run_queue(env, site, submit_fq_count_wdl):
    """ tests that the jaws run queue command has the run id in the stdout."""
    env = env
    site = site

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run queue | grep '\"id\":' | awk '{print $2}' | tr -d ','" % (env)
    (o,e,r) = pf.submit_cmd(cmd)

    ids=o.split()
    run_id = str(submit_fq_count_wdl['run_id'])
    assert run_id in ids


def test_jaws_run_history(env, site, submit_fq_count_wdl):
    """ tests that the jaws run history command has the run id in the stdout."""
    env = env
    site = site
    run_id = str(submit_fq_count_wdl['run_id'])
    wait_for_run(env,run_id)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run history | grep '\"id\": %s' | awk '{print $2}' | tr -d ','" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert o


def test_jaws_wdl_add(env):
    """ tests that the jaws wdl add command added something."""

    wdl= """workflow fq_count { 
         File fastq_file
         call count_seqs { input: infile = fastq_file }
         output { File outfile = count_seqs.outfile }   
        }

        task count_seqs {
        File infile
        command <<< echo ~{infile} >>> 
        output { File outfile = stdout() }
        }
        """

    # write a temporary wdl
    with open(tmp_wdl,"w") as f:
        f.write(wdl)

    # write a temporary readme
    readme='this workflow does not do anything'
    with open(tmp_readme,"w") as f:
        f.write(readme)

    # make sure this wdl doesn't already exist in the catalog (i.e. delete failed for the last test).
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list | grep %s" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    # if command succeeds, then there is still an old wdl in catalog.  We need to delete it.
    if r == 0:
        # delete wdl from catalog so we can test adding it back again
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl delete %s 2.0.0" % (env,wdl_catalog_name)
        (o,e,r) = pf.submit_cmd(cmd)

        assert not r
            
    # add to catalog
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl add %s 2.0.0 %s %s" % (env,wdl_catalog_name,tmp_wdl,tmp_readme)
    (o,e,r) = pf.submit_cmd(cmd)

    # show that it was added
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list | grep %s" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert o

def test_jaws_wdl_update_readme(env):
    """Update a readme for a WDL in the JAWS catalog"""
    readme = 'this readme has been changed'
    with open(tmp_readme,"w") as f:
        f.write(readme)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-doc %s 1.0.0 %s" % (env,released_wdl_catalog_name,tmp_readme)
    (o,e,r) = pf.submit_cmd(cmd)
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl about %s 1.0.0 " % (env,released_wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    assert readme in o

def test_jaws_wdl_update_wdl(env):
    """Update a WDL from the JAWS catalog"""
    readme = 'this readme has been changed'
    with open(tmp_readme,"w") as f:
        f.write(readme)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-doc %s 2.0.0 %s" % (env,wdl_catalog_name,tmp_readme)
    (o,e,r) = pf.submit_cmd(cmd)
    
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl about %s 2.0.0 " % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    assert readme in o

def test_jaws_wdl_update_released_wdl(env):
    """You should not be able to changed a WDL that has been released"""

    task = """task secondecho {
    command { echo second task }
    }"""

    with open(tmp_wdl,"a") as f:
        f.write(task)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-wdl %s 1.0.0 %s" % (env,released_wdl_catalog_name,tmp_wdl)
    (o,e,r) = pf.submit_cmd(cmd)

    assert "Action not allowed" in e

def test_jaws_wdl_versions(env):
    """Test that we can the version for a given WDL
    "fq_count:1.0.0": {
        "created": "2020-11-04T22:24:11Z",
        "last_updated": "2020-11-04T22:24:11Z",
        "name": "fq_count",
        "owner": "jfroula",
        "production_release": "no",
        "version": "1.0.0"
        }
    """
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl versions %s | grep version | awk '{print $2}' " % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    assert o.strip() == '\"2.0.0\"'

def test_jaws_wdl_delete(env):
    """Check that a WDL is deleted"""
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl delete %s 2.0.0" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list | grep %s" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)

    # if grep found wdl_catalog_name still in the catalog, delete failed
    assert not o

def test_jaws_wdl_metadata(env, site, submit_fq_count_wdl):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    env = env
    site = site

    run_id = str(submit_fq_count_wdl['run_id'])
    wait_for_run(env,run_id)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run metadata %s | grep workflowRoot | awk '{print $2}' | tr -d '\"'" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    assert o 

def test_jaws_wdl_errors(env, site, submit_fq_count_wdl):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    env = env
    site = site

    run_id = str(submit_fq_count_wdl['run_id'])
    wait_for_run(env,run_id)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run errors %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    # we don't have any errors in our wdl so the errors command should't return anything.
    # We don't have a good test for this command here, but we only test that the return code is 0.
    assert r == 0
        
def test_jaws_wdl_task_status(env, site, submit_fq_count_wdl):
    """Check that jaws run task-status returns something like this:
     fq_count.count_seqs 1   25177   running success 2021-01-13 12:37:45     The job completed successfully
    """
    env = env
    site = site

    run_id = str(submit_fq_count_wdl['run_id'])
    wait_for_run(env,run_id)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-status %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    assert 'fq_count.count_seqs' in o and 'The job completed successfully' in o

def test_jaws_wdl_log(env, site, submit_fq_count_wdl):
    """Check that the final line of jaws run log returns something like this:
        downloading download complete   2021-01-13 12:41:28 
    """
    env = env
    site = site

    run_id = str(submit_fq_count_wdl['run_id'])
    wait_for_run(env,run_id)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run log %s | tail -1" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    assert 'download complete' in o

def test_jaws_wdl_task_log(env, site, submit_fq_count_wdl):
    """Check that the final line of jaws run log returns something like this:
        downloading download complete   2021-01-13 12:41:28 
    """
    env = env
    site = site

    run_id = str(submit_fq_count_wdl['run_id'])
    wait_for_run(env,run_id)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-log %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    a=[]
    for line in o.split("\n"):
        a.append(line.split())

    # remove empty elements
    a = list(filter(None,a))

    assert a[1][3] == 'created' and a[-1][4] == 'success'

