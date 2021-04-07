#!/usr/bin/env python
"""These functions are to test the "testcases" from the "score_card" integration tests.
google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
There is no need to import conftest.py as it is done automatically.

Author: Angie Kollmer <akollmer@lbl.gov>
Updated: 03/04/21
"""

import os
import json
import pytest
import submission_utils as util

check_tries = 50  # try this many times when waiting for a JAWS run to complete.
check_sleep = 30  # wait for this amount of time between tries.

def test_json_file_does_not_exist(env,dir,site):
    # Submit job with path to json file that does not exist
    # Can't use pf.submit_one_run_to_env here because it exits if submission not successful
    wdl = os.path.join(dir,"WDLs/fq_count.wdl")
    json = "./FileDoesNotExist.json"

    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    submit_cmd = "jaws run submit %s %s %s" % (wdl, json, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)
    #print("cmd: %s\nout: %s\nerror: %s", cmd, o, e)

    # check for the correct error message
    assert "No such file or directory:" in e


def test_json_bad_path_to_input_file_msg(env,dir,site):
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = os.path.join(dir,"WDLs/fq_count.wdl")

    # Submit job with json contains path to a non-existent input file
    # Can't use pf.submit_one_run_to_env here because it exits if submission not successful
    json = "test-inputs/bad_path_inputs.json"
    submit_cmd = "jaws run submit %s %s %s" % (wdl, json, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)
    #print("cmd: %s\nout: %s\nerror: %s" % (cmd, o, e.replace('\n',' ')))

    # check for the correct error message
    assert "UserWarning: Input path not found:" in e

# TODO add more tests to cover the other scorecard scenarios that I said this test would cover
