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
import submission_utils as util

check_tries = 50  # try this many times when waiting for a JAWS run to complete.
check_sleep = 30  # wait for this amount of time between tries.


def test_json_file_does_not_exist(env, dir, site):
    # TESTCASE-4
    # Submit job with path to json file that does not exist
    # Can't use submission_utils submit_wdl function here because it exits if submission not successful
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")
    inputs = os.path.join(dir, "./FileDoesNotExist.json")

    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    submit_cmd = "jaws submit --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "No such file or directory:" in e


def test_input_file_is_not_json_format(env, dir, site):
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")
    # testing for message when inputs file is not json, so using wdl again instead of a json file
    inputs = wdl
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    submit_cmd = "jaws submit --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "is not a valid JSON file" in e


def test_json_bad_path_to_input_file_msg(env, dir, site):
    # TESTCASE-5a
    # Submit job with json that contains a path to a non-existent input file
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")

    # Can't use  submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/bad_path_inputs.json")
    submit_cmd = "jaws submit --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "Input path not found or inaccessible:" in e


def test_misspelled_variable_in_input_file_msg(env, dir, site):
    # TESTCASE-5b
    # Submit job with json that contains a misspelled variable name
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")
    input_json = os.path.join(dir, "test-inputs/misspelled_variable.json")

    # we CAN use submission utils here because this job submits successfully
    # error isn't seen until the run fails
    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    run_id = data['run_id']
    util.wait_for_run(run_id, env, check_tries, check_sleep)

    # check for the correct error message
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    errors_cmd = "jaws errors %s" % (run_id)
    cmd = source_cmd + errors_cmd
    (r, o, e) = util.run(cmd)

    # in 2.2 this error is seen in the metadata
    # I think that in 2.3 the error should also be displayed by the error command output
    assert "Required workflow input 'fq_count.fastq_file' not specified" in o


def test_bad_input_file_permissions_msg(env, dir, site):
    # TESTCASE-6
    # Submit json that contains a path to a file with bad permissions
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = os.path.join(dir, "WDLs/fq_count.wdl")

    # Can't use submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/bad_permissions.json")
    submit_cmd = "jaws submit --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "Input path not found or inaccessible:" in e


def test_invalid_wdl_syntax_msg(env, dir, site):
    # TESTCASE-7
    # Submit invalid WDL syntax
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = os.path.join(dir, "WDLs/bad_syntax.wdl")

    # Can't use submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/fq_count.json")
    submit_cmd = "jaws submit --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "ERROR: Unexpected symbol" in e


def test_invalid_wdl_semantics_msg(env, dir, site):
    # TESTCASE-8
    # Submit invalid WDL semantics
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = os.path.join(dir, "WDLs/bad_semantics.wdl")

    # Can't use submission_utils submit_wdl here because it exits if submission not successful
    inputs = os.path.join(dir, "test-inputs/fq_count.json")
    submit_cmd = "jaws submit --no-cache %s %s %s" % (wdl, inputs, site)
    cmd = source_cmd + submit_cmd
    (r, o, e) = util.run(cmd)

    # check for the correct error message
    assert "No input xxx_file found" in e
