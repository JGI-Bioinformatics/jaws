#!/usr/bin/env ipython
"""
These functions are to test the "testcases" from the "score_card" integration tests.
google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
"""

import os
import pytest
import json
import time
import submission_utils as util


WDL = "/WDLs/timeout.wdl"
INP = "/test-inputs/timeout.json"

check_sleep = 30
check_tries = 50


def test_timeout(env, dir, site):
    wdl = dir + WDL
    input_json = dir + INP

    run_id = util.submit_wdl(env, wdl, input_json, site)["run_id"]
    util.wait_for_run(run_id, env, check_tries, check_sleep)

    time.sleep(60)

    # get the errors from JAWS for that run
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    errors_cmd = "jaws errors %s" % (run_id)
    cmd = source_cmd + errors_cmd
    r, o, e = util.run(cmd) 

    # do the check!
    fail_msg = "error. Keyword absent: \"timeout\" (%s)" % run_id
    assert "failed with timeout" in o, fail_msg


def test_scatter_timeout(env, submit_scatter_timeout):
    """
    TESTCASE-44
    When user submits a wdl with a scatter function and a timeout occurs in the scatter function
    then the timeout message should appear in the output from the errors command
    """
    run_id = str(submit_scatter_timeout["run_id"])
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws errors %s" % (env, run_id)
    (r, o, e) = util.run(cmd)

    assert "failed with timeout" in o, "scatter timeout error should appear in errors"
