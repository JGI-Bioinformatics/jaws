#!/usr/bin/env python

# These functions are to test the "score_card" unit tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451
#
# Specifically, this script will make sure job submissions in JAWS returns run 
# statistics that a user might desire, like which WDL task is being run.
# Tests covered are TESTCASE2 & 3
#
# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
#
# Test Cases as defined in the above google doc:
#  TESTCASE-2
#  run level info
#  With a Run-ID the user should be able to receive the following information through the JAWS CLI about the run during the run:
#  * State (eg ready, uploading, submitted, running, succeeded, downloading, finished, failed) => jaws run status
#  * Information since when run is in said state jaws run status: updated
#  * A list of transitions between states (entered state at time, left state at time) =>  jaws run log
#  * A list of tasks with unique IDs for the run => jaws run task-status & task-log
#  * Output data for the run (shows paths) => jaws run status
#  * Input data for the run (shows paths) -- inputs json copied to output directory at beginning of run
#  * WDL specification used for the run -- wdl copied to output directory at beginning of run
#
#  TESTCASE-3
#  Task level info
#  With a Run-ID the user should be able to receive the following information through the JAWS CLI about the task:
#  * has state info jaws run task-status
#  * timestamp since when task was in state jaws run task-status
#  * A list of transitions between states (entered state at time, left state at time)  jaws run task-log
#  * stderr/stdout of the task output


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
def test_jaws_run_log(env,submit_wdl_and_wait):
    """ * State (eg ready, uploading, submitted, running, succeeded, downloading, finished, failed) => jaws run log
        * Information since when run is in said state jaws run status: updated
        * A list of transitions between states (entered state at time, left state at time) =>  jaws run log
        * Output data for the run (shows paths) => jaws run status

        results from log
        ------------------
        #STATUS_FROM	STATUS_TO	TIMESTAMP	REASON
        created	        uploading	   2021-02-02 22:06:10	upload_task_id=dba58112-65a2-11eb-827f-0275e0cda761
        uploading	upload complete	   2021-02-02 22:06:21
        upload          complete	   2021-02-02 22:06:33	cromwell_run_id=2c838624-0670-4446-a0ec-9caa451e723f
        submitted	queued	           2021-02-02 22:06:45
        queued	        running	           2021-02-02 22:07:39
        running	        succeeded	   2021-02-02 22:07:55
        succeeded	downloading	   2021-02-02 22:08:07	download_task_id=20ceb2ea-65a3-11eb-8c38-0eb1aa8d4337
        downloading	download complete  2021-02-02 22:08:54
    """

    run_id = submit_wdl_and_wait['run_id']
    
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run log %s | tail -n+2" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    stages = []
    line_list = o.split("\n")
    line_list = list(filter(None, line_list))  # remove empty element
    for i in line_list:
        stages.append(i.split("\t")[1]) # index 1 should be the STATUS_TO column

    # do we have a log of all the states
    expected=['uploading', 'upload complete', 'submitted', 'queued', 'running', 'succeeded', 'downloading', 'download complete']
    a=0
    for stage in expected:
        if stage in stages:
            a += 1

    if a == 8:
        assert 1
    else:
        assert 0


    # do we have timestamps
    for i in line_list:
        times = i.split("\t")[2].split()

        # We should have this format: 2021-02-02 22:07:55
        if re.search(r"^\d{4}-\d{2}-\d{2}$", times[0]):
            assert 1
        else:
            assert 0

        if re.search(r"^\d{2}:\d{2}:\d{2}$", times[1]):
            assert 1
        else:
            assert 0

    # test that an output directory was displayed 
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    data = json.loads(o)
    assert os.path.exists(data['output_dir'])
    
def test_jaws_run_task_log(env,submit_wdl_and_wait):
    """ 
        TESTCASE-3
        Task level info
        With a Run-ID the user should be able to receive the following information through the JAWS CLI about the task:
        * has state info jaws run task-status
        * timestamp since when task was in state jaws run task-status
        * A list of transitions between states (entered state at time, left state at time)  jaws run task-log

        results from task-log
        ---------------------
        #TASK_NAME	ATTEMPT	CROMWELL_JOB_ID	STATUS_FROM	STATUS_TO	TIMESTAMP	REASON
        fq_count.count_seqs	1	44472	created	ready	2021-02-02 22:06:39
        fq_count.count_seqs	1	44472	ready	queued	2021-02-02 22:06:41
        fq_count.count_seqs	1	44472	queued	pending	2021-02-02 22:06:42
        fq_count.count_seqs	1	44472	pending	running	2021-02-02 22:07:32
        fq_count.count_seqs	1	44472	running	success	2021-02-02 22:07:34
    """


    run_id = submit_wdl_and_wait['run_id']
    
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-log %s | tail -n+2" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    stages_from = []
    stages_to = []
    times = []

    line_list = o.split("\n")
    line_list = list(filter(None, line_list))  # remove empty element
    for i in line_list:
        stages_from.append(i.split("\t")[3]) # index 3 should be the STATUS_FROM column
        stages_to.append(i.split("\t")[4]) # index 4 should be the STATUS_TO column
        times.append(i.split("\t")[5]) # index 5 should be the TIMESTAMP column

    # do we have a log of all the states
    status_from_expected=['created', 'ready', 'queued', 'pending', 'running']
    status_to_expected=['ready', 'queued', 'pending', 'running', 'success']

    for stage in status_from_expected:
        if stage not in stages_from:
            assert 0

    for stage in status_to_expected:
        if stage not in stages_to:
            assert 0

    # do we have timestamps
    for date_and_time in times:
        date = date_and_time.split()[0]
        time = date_and_time.split()[1]

        # We should have this format: 2021-02-02 22:07:55
        if re.search(r"^\d{4}-\d{2}-\d{2}$", date):
            assert 1
        else:
            assert 0

        if re.search(r"^\d{2}:\d{2}:\d{2}$", time):
            assert 1
        else:
            assert 0

