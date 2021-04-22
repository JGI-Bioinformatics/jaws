#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory. There is no need to import conftest.py as it is done automatically.

import sys
import os
import re
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
#
#
def test_jaws_run_task_log(env,submit_fq_count_wdl):
    """ 
    1) I will check that the input WDL and json file are saved in the output dir.

    2) I will also check that the raw cromwell file structure was created

    3) I will check that the expected files were created in the execution dir (i.e. stdout and stderr)
        This cromwell file structure should exist
        fq_count_out/call-count_seqs/execution/
            num_seqs.txt  rc  script  script.submit  stderr  stderr.submit	stdout	stdout.submit
    """
    run_id = str(submit_fq_count_wdl['run_id'])
    output_dir = submit_fq_count_wdl['output_dir']

    util.wait_for_run(run_id,env,check_tries,check_sleep)

    # grab submission_id from "jaws run status" command so we can construct name of output files
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env,run_id)
    (r,o,e) = util.run(cmd)
    status_data = json.loads(o)
    submission_id = status_data['submission_id']
    input_wdl = submission_id + ".wdl"
    input_json = submission_id + ".json"

    # check that we have the initial WDL saved to the output_dir
    # using the full path (output_dir and input_wdl), we are essentially testing that the output_dir 
    # was correct and that the wdl file got created.

    # verify it is a valid wdl
    with open(os.path.join(output_dir,input_wdl)) as fh:
        if not "workflow fq_count" in fh.readline():
            assert 0,"This does not look like a valid workflow"

    # check that we have a valid inputs json
    with open(os.path.join(output_dir,input_json)) as fh:
        expected = '"fq_count.fastq_file":' 
        if not expected in fh.read():
            assert 0,"This does not look like a valid inputs json file"
    
    expected_files = ["num_seqs.txt","rc","script","script.submit","stderr","stderr.submit","stdout","stdout.submit"]
    for file in expected_files:
        if not os.path.exists(os.path.join(output_dir,"call-count_seqs/execution/",file)):
            assert 0, (f"expected result file not found in: {os.path.join(output_dir,'call-count_seqs/execution/',file)}")

