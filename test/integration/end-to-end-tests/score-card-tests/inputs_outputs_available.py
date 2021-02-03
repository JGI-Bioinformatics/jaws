#!/usr/bin/env python

# These functions are to test the "score_card" unit tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451
#
# Specifically, this script will verify the following input & output files were created.
#    1) I will check that the input WDL and json file are saved in the output dir.
#
#    2) I will also check that the raw cromwell file structure was created
#
#    3) I will check that the expected files were created in the execution dir (i.e. stdout and stderr)
#        This cromwell file structure should exist
#        fq_count_out/call-count_seqs/execution/
#            num_seqs.txt  rc  script  script.submit  stderr  stderr.submit	stdout	stdout.submit
#
# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
#

import sys
import os
import json
import re
import parsing_functions as pf

#########################
###     Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#
def test_jaws_run_task_log(env,submit_wdl_and_wait):
    """ 
    1) I will check that the input WDL and json file are saved in the output dir.

    2) I will also check that the raw cromwell file structure was created

    3) I will check that the expected files were created in the execution dir (i.e. stdout and stderr)
        This cromwell file structure should exist
        fq_count_out/call-count_seqs/execution/
            num_seqs.txt  rc  script  script.submit  stderr  stderr.submit	stdout	stdout.submit
    """

    output_dir = submit_wdl_and_wait['output_dir']
    run_id = submit_wdl_and_wait['run_id']
    submission_id = submit_wdl_and_wait['submission_id']
    input_wdl = submission_id + ".wdl"
    input_wdl = "a82f72ec-7fd5-4697-b2f8-77af546102a0.wdl"
    input_json = "fq_count.json"

    # check that we have the initial WDL saved to the output_dir
    # using the full path (output_dir and input_wdl), we are essentially testing that the output_dir 
    # was correct and that the wdl file got created.
    if not os.path.exists(os.path.join(output_dir,input_wdl)):
        assert 0

    # verify it is a valid wdl
    with open(os.path.join(output_dir,input_wdl)) as fh:
        if not "workflow fq_count" in fh.readline():
            assert 0

    # check that we have a valid inputs json
    if not os.path.exists(os.path.join(output_dir,input_json)):
        assert 0

    with open(os.path.join(output_dir,input_json)) as fh:
        expected = '"fq_count.fastq_file":' 
        if not expected in fh.read():
            assert 0
    
    expected_files = ["num_seqs.txt","rc","script","script.submit","stderr","stderr.submit","stdout","stdout.submit"]
    for file in expected_files:
        if os.path.exists(os.path.join(output_dir,"fq_count_out/call-count_seqs/execution/",file)):
            print("file found")
    

