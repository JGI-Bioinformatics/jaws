#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory. There is no need to import conftest.py as it is done automatically.

import pytest
import json
import time
import submission_utils as util

# set variables specific for this series of tests
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
# This test requires "env" which is from the command line and captured by a fixture inside the conftest.py file.
# For example --env prod on the command line passes "prod" to this env variable.
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

# This test requires a fixture "submit_fq_count_wdl" which is inside the conftest.py file.
# This variable contains the json stdout that is created when you do a WDL JAWS submission.
# Also, " --env prod" on the command line passes "prod" to this env variable.
def test_jaws_run_history(env, submit_fq_count_wdl):
    """ tests that the jaws run history command has the run id in the stdout."""
    data = submit_fq_count_wdl
    run_id = str(data['run_id'])
    util.wait_for_run(env,run_id,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run history | grep '\"id\": %s' | awk '{print $2}' | tr -d ','" % (env,run_id)
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert o


