#!/usr/bin/env python

# These functions are to test the "score_card" unit tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451
#
# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
#
# The following tests will make sure JAWS handles subworkflows correctly.
# 1) task-status verifies all subworkflows task status was shown
# 2) raw cromwell subworkflow files are returned to user defined output dir
# 3) subworkflow WDLs are saved in the user defined output dir
# 5) metadata command also returns cromwell metadata for subworkflows

import sys
import os
import json
import re
import pytest
import time
import configparser
import glob
import parsing_functions as pf


# get some environmental vars
config = configparser.ConfigParser()
config.read(os.environ.get('MYINI_FILE'))

WDL        = config['wdl']['wdl']
INPUT_JSON = config['wdl']['input_json']
ENV        = config['wdl']['env']
OUTDIR     = config['wdl']['outdir']
SITE       = config['wdl']['site']

#########################
###     Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#
def test_task_log(submit_wdl_and_wait):
    """
    # task-status verifies all subworkflows task status was shown
    #
    #TASK_NAME	ATTEMPT	CROMWELL_JOB_ID	STATUS_FROM	STATUS_TO	TIMESTAMP	REASON	STATUS_DETAIL
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard      1	46806	running	success	2021-02-08 20:53:55  The job completed successfully
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing	1	46807	running	success	2021-02-08 20:53:58  The job completed successfully
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment  1	46808	running	success	2021-02-08 20:54:25  The job completed successfully
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams 1	46809	running	success	2021-02-08 20:54:37  The job completed successfully
    main_wdl.bam_stats	                              1	46810	running	success	2021-02-08 20:56:36  The job completed successfully
    """
    run_id = submit_wdl_and_wait['run_id']
    time.sleep(30)  # wait some time before running task-status since there is some lag between 
                    # when "jaws run status" calls success and when "jaws run task-status" calls success.
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-status %s | tail -n+2" % (ENV,run_id)
    (o,e,r) = pf.submit_cmd(cmd)

    # put the table into a dictionary
    task_names = []
    status_to = []
    line_list = o.split("\n")
    line_list = list(filter(None, line_list))  # remove empty element
    for i in line_list:
        task_names.append(i.split("\t")[0])
        status_to.append(i.split("\t")[4])

    # check that the subworkflows tasks are in the list
    assert len(task_names) == 5

    # make sure all tasks completed with success
    assert len(list(filter(lambda x: (x == 'success'),status_to))) == 5
    
def test_for_raw_cromwell_files():
    """ 
    test that raw cromwell subworkflow files are returned to user defined output dir.

    These files should exist OUTDIR/call-bbmap_shard_wf/<cromwell-hash>/
        shard
        bbmap_indexing
        alignment
        merge_bams

        out/call-bbmap_shard_wf/align.bbmap_shard_wf/c9f67f71-1acd-4d8c-8879-14b7b1a28b54/call-shard/execution/rc
    """
    cmd = "find %s/call-bbmap_shard_wf/align.bbmap_shard_wf -name rc -exec cat {} \\; | grep -c 0" % (OUTDIR)
    (o,e,r) = pf.submit_cmd(cmd)

    # make sure all 4 of our "rc" files returned 0
    assert int(o.strip()) == 4


def test_saved_subwdl(submit_wdl_and_wait):
    """
    subworkflow WDLs are saved in the user defined output dir

    """
    submission_id = submit_wdl_and_wait['submission_id']
    zip_file = os.path.join(OUTDIR,submission_id + ".zip")

    assert os.path.exists(zip_file)

    cmd = "unzip -l %s" % (zip_file)
    (o,e,r) = pf.submit_cmd(cmd)
    assert "alignment.wdl" in o


def test_subworkflow_metadata(submit_wdl_and_wait):
    """
    metadata command also returns cromwell metadata for subworkflows
    """
    run_id = submit_wdl_and_wait['run_id']
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run metadata %s" % (ENV,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    meta_output = json.loads(o)

    # make sure metadata returns these calls
    expected = ['bbmap_shard_wf.alignment', 'bbmap_shard_wf.bbmap_indexing', 'bbmap_shard_wf.merge_bams', 'bbmap_shard_wf.shard']
    calls=[]
    for key in meta_output:
        for yek in meta_output[key]['calls']:
            calls.append(yek)
        
    assert len([x for x in expected if x in calls]) == 4
        

    

