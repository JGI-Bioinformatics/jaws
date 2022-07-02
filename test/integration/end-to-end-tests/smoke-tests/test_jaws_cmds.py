#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n  # noqa
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.

import os
import json
import uuid
import shutil
import submission_utils as util
import pytest

check_tries = 360
check_sleep = 60


#####################
#     Functions     #
#####################
#
# Test functions for verification of jaws log commands (log,task-log,status).
#
def test_jaws_info():
    """tests that there is a valid output for jaws info. Name should be dev,staging, or prod and version should have some value.
    {
    "docs_url": "https://jaws-docs.readthedocs.io/en/latest/",
        "name": "prod",
        "version": "2.1"
    }"""  # noqa

    cmd = "jaws info"
    (r, o, e) = util.run(cmd)

    data = json.loads(o)

    # do we have an acceptable name
    assert data["name"] in ["prod", "staging", "dev"]
    assert data["version"] is not None


def test_jaws_health():
    """tests that the jaws health is working. We don't care if some services are down.
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

    cmd = "jaws health"
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    actual_keys = list(data.keys())
    required_keys = [
        "CORI-Cromwell",
        "CORI-RMQ",
        "CORI-Site",
        "JAWS-Central",
        "JGI-Cromwell",
        "JGI-RMQ",
        "JGI-Site",
    ]

    for k in required_keys:
        assert k in actual_keys

def test_jaws_status(submit_fq_count_wdl):
    """test that the jaws status --verbose command has some of the expected keys displayed in its stdout"""
    data = submit_fq_count_wdl
    run_id = str(data["run_id"])
    cmd = "jaws status --verbose %s" % (run_id)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)
    assert "compute_site_id" in data
    assert "id" in data
    assert "input_site_id" in data
    assert "json_file" in data
    assert "result" in data
    assert "status" in data
    assert "submitted" in data
    assert "tag" in data
    assert "updated" in data
    assert "user_id" in data

def test_jaws_queue(site, dir):
    """tests that the jaws queue command has the correct site in the output when --site is used."""

    wdl = dir + "/WDLs/fq_count.wdl"
    myjson = dir + "/test-inputs/fq_count.json"
    cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, myjson, site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)
    run_id = data["run_id"]

    # jaws queue --site [CORI, JGI, TAHOMA]
    cmd = "jaws queue --site %s" % (site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    # if the site_id is found in the jaws queue output, I don't want to
    # assert true yet, I want to cancel the run and then I can use assert. So
    # just save the boolean in 'result' for now.
    result = False
    has_id = False
    ids = []
    if data:
        for d in data:
            if d["compute_site_id"].lower() == site.lower():
                result = True
                ids.append(d["id"])
    else:
        assert result, f"no runs were found in the queue for site: {site}"

    if run_id in ids:
        has_id = True

    cmd = "jaws cancel --quiet %s" % (run_id)
    (r, o, e) = util.run(cmd)

    assert result and has_id


def test_jaws_history(submit_fq_count_wdl):
    """tests that the jaws history command has the run id in the stdout."""
    data = submit_fq_count_wdl
    run_id = str(data["run_id"])
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws history | grep '\"id\": %s' | awk '{print $2}' | tr -d ','" % (run_id)
    (r, o, e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert o


def test_jaws_metadata(site, submit_fq_count_wdl):
    """Check that a jaws metadata returns workflowRoot has a value"""

    # skip this fixture if run on aws
    if "aws" in site.lower():
        pytest.skip("this test won't work on AWS")

    data = submit_fq_count_wdl
    run_id = str(data["run_id"])
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws metadata %s | grep workflowRoot | awk '{print $2}' | tr -d '\"'" % (
        run_id
    )
    (r, o, e) = util.run(cmd)

    assert o

@pytest.mark.xfail(reason="we created ticket for bug in errors command")
def test_jaws_errors(submit_bad_task):
    """Check that jaws errors catches the stderr error"""
    run_id = str(submit_bad_task["run_id"])
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws errors %s" % (run_id)
    (r, o, e) = util.run(cmd)

    assert "command not found" in o


def test_jaws_log(submit_fq_count_wdl):
    """Check that the first line of jaws log returns something like this:
    #STATUS_FROM       STATUS_TO          TIMESTAMP            COMMENT
    created            upload queued      2022-05-24 16:32:41
    upload queued      upload complete    2022-05-24 16:32:52
    upload complete    ready              2022-05-24 16:33:03
    ready              submitted          2022-05-24 16:33:45  cromwell_run_id=1c222955-b4c2-4315-b3fd-7c4f20ffaa2d
    submitted          queued             2022-05-24 16:33:56
    queued             running            2022-05-24 16:37:23
    running            succeeded          2022-05-24 16:37:23
    succeeded          finished           2022-05-24 16:37:34
    finished           download queued    2022-05-24 16:37:40
    download queued    download complete  2022-05-24 16:39:16
    download complete  email sent         2022-05-24 16:39:28
    email sent         done               2022-05-24 16:39:39
    """
    run_id = str(submit_fq_count_wdl["run_id"])
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws log %s " % (run_id)
    (r, o, e) = util.run(cmd)
    a = o.split("\n")
    a = list(filter(None, a))  # remove empty elements

    assert "created" in a[1] and "upload queued" in a[1]

    # there should be 4 columns of output (but time is split into two so actually 5)
    assert len(a[1].split()) == 5

    # there should be 13 output rows (including header)
    assert len(a) == 13 or len(a) == 11 ## We added 11 because we are missing the last step "email sent --> done". V2.8.6


def test_jaws_get(submit_bad_task):
    """
    Check that the 'get' cmd works, even when a run fails.
    """
    run_id = str(submit_bad_task["run_id"])

    # remove any old copy of the temp dir if exists
    outdir = str(uuid.uuid4())

    cmd = "jaws get --quiet --complete %s %s" % (run_id, outdir)
    (r, o, e) = util.run(cmd)
    assert r == 0

    # These should exist: (script, script.submit, stderr.submit, stdout.submit)
    assert os.path.exists(
        os.path.join(outdir, "call-count_seqs/execution/script")
    )

    try:
        shutil.rmtree(outdir)
    except OSError as error:
        print(f"Error: {outdir}: {error}")


def test_tag(submit_fq_count_wdl):
    """
    Check that the '--tag' flag is showing 'submit_fq_count_wdl'.
    """
    run_id = str(submit_fq_count_wdl["run_id"])
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws status %s" % (run_id)
    (r, o, e) = util.run(cmd)
    assert r == 0
    data = json.loads(o)

    assert data["tag"] == "submit_fq_count_wdl"


def test_jaws_history_site_filter(site, submit_fq_count_wdl):
    """
    jaws history --site [CORI, JGI, TAHOMA]
    """
    site = "cori"
    cmd = "jaws history --site %s" % (site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["compute_site_id"].lower() == site.lower()
    else:
        assert 0, f"no runs were found in the history for site: {site}"


def test_jaws_history_result_filter_succeeded(submit_fq_count_wdl):
    """
    jaws history --result [succeeded, failed]
    Checking the output only with "succeeded" and "failed"
    """

    run_id = submit_fq_count_wdl["run_id"]
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws history --result succeeded"
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["result"] == "succeeded"
    else:
        assert 1, "no runs were found in the history that have result: succeeded"


def test_jaws_history_result_filter_failed(submit_bad_task):
    """
    jaws history --result failed
    Checking the output only with "succeeded" and "failed"
    """
    run_id = submit_bad_task["run_id"]
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws history --result failed"
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["result"] == "failed"
    else:
        assert 1, "no runs were found in the history that have result: failed"


#def test_jaws_task_summary(submit_fq_count_wdl):
#    """
#    jaws task-summary <run_id>
#
#    #NAME                CROMWELL_JOB_ID  CACHED  RESULT   QUEUED               QUEUE_WAIT  RUNTIME  MAX_TIME
#    fq_count.count_seqs  158914           False   success  2022-01-24 20:50:21  0:12:44     0:00:01  00:10:00
#
#    Check that the fourth column is "success" 
#    """
#    run_id = str(submit_fq_count_wdl["run_id"])
#    util.wait_for_run(run_id, check_tries, check_sleep)
#
#    cmd = (f"jaws task-summary --fmt json {run_id}")
#    (r, o, e) = util.run(cmd)
#    data = json.loads(o)
#
#    if data:
#        assert data[0][3] == 'success'
#    else:
#        assert 1, "task-summary did not return \"success\" in the fourth column" 
