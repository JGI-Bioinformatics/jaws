#!/usr/bin/env python
"""
These functions are to test the "testcases" from the "score_card" integration tests.
google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

This library of tests uses "fixtures" from conftest.py which should be located in the same directory. 
There is no need to import conftest.py as it is done automatically.

This script will test that fq_count.wdl & fq_count.json will successfully complete on JAWS on [cori|jgi] and [staging|prod].
It is to be run as a nightly health monitoring test.

Author: Jeff Froula <jlfroula@lbl.gov>
Created: 03/24/21
"""

import pytest
import json
import time
import submission_utils as util

# set variables specific for this series of tests
check_tries=100
check_sleep=30

#########################
###     Functions     ###
#########################
def test_success_at_staging_cori():
    """ the submit_and_wait_for_success_prod.v3.2.sh command is for testing staging on cori and jgi. It will submit fq_count.wdl and send the result to slack: jaws_health """

    cmd = "./submit_and_wait_for_success_prod.v3.2.sh -r staging -s cori"
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert "UP" in o 

def mtest_success_at_staging_jgi():
    """ the submit_and_wait_for_success_prod.v3.2.sh command is for testing staging on cori and jgi. It will submit fq_count.wdl and send the result to slack: jaws_health """

    cmd = "./submit_and_wait_for_success_prod.v3.2.sh -r staging -s jgi"
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert "UP" in o 

def mtest_success_at_prod_cori():
    """ the submit_and_wait_for_success_prod.v3.sh command is for testing prod on cori and jgi. It will submit fq_count.wdl and send the result to slack: jaws_health """

    cmd = "./submit_and_wait_for_success_prod.v3.sh -r prod -s cori"
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert "UP" in o

def mtest_success_at_prod_jgi():
    """ the submit_and_wait_for_success_prod.v3.sh command is for testing prod on cori and jgi. It will submit fq_count.wdl and send the result to slack: jaws_health """

    cmd = "./submit_and_wait_for_success_prod.v3.sh -r prod -s jgi"
    (r,o,e) = util.run(cmd)

    # if there was something in stdout, then grep found "id": <run_id>
    assert "UP" in o

