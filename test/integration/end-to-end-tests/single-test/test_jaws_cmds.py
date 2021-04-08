#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.

import os,sys
import pytest
import json
import time
import submission_utils as util

tmp_wdl = "pow23.wdl"
tmp_readme = "pow23.md"
wdl_catalog_name="tmp_wdl_catalog_name"
released_wdl_catalog_name="fq_count"
check_tries=100
check_sleep=30


#########################
###     Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#
def test_jaws_info():
    """ tests that there is a valid output for jaws info. Name should be dev,staging, or prod and version should have some value.
    {
    "docs_url": "https://jaws-docs.readthedocs.io/en/latest/",
        "name": "prod",
        "version": "2.1"
    } """

    cmd = "jaws info"
    (r,o,e) = util.run(cmd)

    data = json.loads(o)

    # do we have an acceptable name
    assert data["name"] in ["prod","staging","dev"]
    assert data["version"] is not None

def test_jaws_status():
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

    cmd = "jaws status"
    (r,o,e) = util.run(cmd)
    data = json.loads(o)

    actual_keys = list(data.keys())
    required_keys = ["CORI-Cromwell","CORI-RMQ","CORI-Site","JAWS-Central","JGI-Cromwell","JGI-RMQ","JGI-Site"]

    for k in required_keys:
        assert k in actual_keys


def test_jaws_run_queue(env,submit_fq_count_wdl):
    """ tests that the jaws run queue command has the run id in the stdout."""

    data = submit_fq_count_wdl
    run_id = str(data['run_id'])

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run queue | grep '\"id\":' | awk '{print $2}' | tr -d ','" % (env)
    (r,o,e) = util.run(cmd)
    ids=o.split()
    assert run_id in ids
