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

check_tries = 360  # try this many times when waiting for a JAWS run to complete.
check_sleep = 60  # wait for this amount of time between tries.


def test_json_file_does_not_exist(dir, site):
    # TESTCASE-4
    # Submit job with path to json file that does not exist
    # Can't use submission_utils submit_wdl function here because it exits if submission not successful
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")
    inputs = os.path.join(dir, "./FileDoesNotExist.json")

    cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, inputs, site)
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "No such file or directory:" in e


def test_input_file_is_not_json_format(dir, site):
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")
    # testing for message when inputs file is not json, so using wdl again instead of a json file
    inputs = wdl
    cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, inputs, site)
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "is not a valid JSON file" in e


def test_json_bad_path_to_input_file_msg(dir, site):
    # TESTCASE-5a
    # Submit job with json that contains a path to a non-existent input file
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")

    # Can't use  submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/bad_path_inputs.json")
    submit_cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = submit_cmd
    (r, o, e) = util.run(cmd)
    data = json.loads(o)
    run_id = data["run_id"]

    # check for the correct error message
    assert "WARNING: Input path not found or inaccessible" in e

    cmd = "jaws cancel %s" % (run_id)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)
    assert 'cancelled' in o


def test_misspelled_variable_in_input_file_msg(dir, site):
    # TESTCASE-5b
    # Submit job with json that contains a misspelled variable name
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")
    input_json = os.path.join(dir, "test-inputs/misspelled_variable.json")

    # we CAN use submission utils here because this job submits successfully
    # error isn't seen until the run fails
    submit_cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, input_json, site)
    (r, o, e) = util.run(submit_cmd)

    assert "KeyError: 'fq_count.fastq_file_misspelled'" in e


def test_bad_input_file_permissions_msg(dir, site):
    # TESTCASE-6
    # Submit json that contains a path to a file with bad permissions
    # This test uses the bad_permissions.json which points to a file that has no read permissions
    # set for owner, group or user
    # JAWS should show user an error message that explains the problem to the user

    wdl = os.path.join(dir, "WDLs/fq_count.wdl")

    # Can't use submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/bad_permissions.json")
    submit_cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, inputs, site)
    (r, o, e) = util.run(submit_cmd)

    # check for the correct error message
    assert "no_read_perms.fastq" in e, "file name should be in error message"
    assert (
        "Permission denied" in e
    ), "permissions problem should be explained in error message"


def test_invalid_wdl_syntax_msg(dir, site):
    # TESTCASE-7
    # Submit invalid WDL syntax
    wdl = os.path.join(dir, "WDLs/bad_syntax.wdl")

    # Can't use submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/fq_count.json")
    cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, inputs, site)
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "ERROR: Unexpected symbol" in e


def test_invalid_wdl_semantics_msg(dir, site):
    # TESTCASE-8
    # Submit invalid WDL semantics
    wdl = os.path.join(dir, "WDLs/bad_semantics.wdl")

    # Can't use submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/fq_count.json")
    cmd = "jaws submit --quiet --no-cache %s %s %s" % (wdl, inputs, site)
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "No input xxx_file found" in e
