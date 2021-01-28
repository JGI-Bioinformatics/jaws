#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
# Sample output from the submit wdl fixture:
#    """data = {'output_dir': '/global/cscratch1/sd/jaws/jfroula/reference_db_out', 
#        'output_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 
#        'run_id': 33, 
#        'site_id': 'NERSC', 
#        'status': 'uploading', 
#        'submission_id': '9134e365-bef7-41de-8073-b9633f9233c0', 
#        'upload_task_id': '608ec02e-b026-11ea-beea-0e716405a293'}
#        """

import sys
import os
import json
from time import sleep
from subprocess import Popen, PIPE

def run(cmd):
    output = Popen(cmd, stdout=PIPE,
             stderr=PIPE, shell=True,
             universal_newlines=True)
             #cwd="/global/cscratch1/sd/jaws/jfroula", universal_newlines=True)

    stdout,stderr=output.communicate()
    rc=output.returncode

    return rc,stdout,stderr


def test_services_are_up():

    cmd = 'jaws status'
    rc,stdout,stderr = run(cmd)
    data = json.loads(stdout)
    assert data['JAWS-Central']   == 'UP'
    assert data['NERSC-Cromwell'] == 'UP'
    assert data['NERSC-RMQ']      == 'UP'
    assert data['NERSC-Site']     == 'UP'

def mytest_status(submit_wdl):
#def test_status():
    """
    This test expects that the wdl has finished processing.
    We test that some basic information is availabe from the "jaws run status" command.

    "cromwell_run_id": "189726a8-b5c6-49de-a1e3-032faae78e3e",
    "download_task_id": "e57ab266-b026-11ea-beea-0e716405a293",
    "id": 33,
    "input_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    "input_site_id": "NERSC",
    "output_dir": "/global/cscratch1/sd/jaws/jfroula/reference_db_out",
    "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    "site_id": "NERSC",
    "status": "finished",
    "status_detail": "The run output has been returned to the user and tmpfiles were purged.",
    "submission_id": "f0875e8b-100c-4c02-9ef8-8606356821d3",
    "submitted": "2020-06-16T23:08:29Z",
    "updated": "2020-06-16 23:13:41",
    "upload_task_id": "4b2e7378-b026-11ea-9a3b-0255d23c44ef"
    """

    # test that we see all expected run stats. The job status will be checked every 10 seconds to 
    # try and capture most of the transition states, but its not practical to expect all will be captured.
    # testcase-2a
    status=''
    list_of_possibles=["uploading","queued", "submitted","running","succeeded","downloading","finished"]
    run_id = str(submit_wdl['run_id'])
    #run_id = "102"
    cmd = 'jaws run status ' + run_id
    status_counts={}

    # Make sure we have most of the expected transition states for a job by 
    # keep checking status of submitted job then wait 20 seconds. 
    end_if_too_long = 0
    while (status != "succeeded" or status != "failed"):
        end_if_too_long += 1
        rc,stdout,stderr = run(cmd)
        data = json.loads(stdout)
        status = data['status']
        print(f"now testing status {status}")
        if status not in list_of_possibles:
            sys.stderr.write(f"status {status} not in acceptable list")
            assert 0
        status_counts={status: list_of_possibles.count(status)}
        print(f"in state {status}")
        sleep(20)
        if end_if_too_long > 10:
            break

    print(f"status counts {status_counts}")
    total=0
    for key,value in status_counts.items():
        total += value

    # test that we saw at least this many states
    assert total > 2

#def test_run_level_info(submit_wdl):
#    output_dir = submit_wdl('output_dir')
#
#    # testcase-2b
##    # test if there is a timestamp for when state began
##    updated = data['updated']
##    s = updated.split()
##    date = s[0].split('-')
##    time = s[1].split(':')
##    assert len(date) == 3
##    assert len(time) == 3
##
#    # testcase-2c
##    # test if there are transition states for each task 
##    # null    ready
##    # ready   queued
##    # queued  pending
##    # pending running
##    # running success
##    cmd = 'jaws run task-log ' + run_id
##    rc,stdout,stderr = run(cmd)
##    data = json.loads(stdout)
##
##    transitions=[]
##    # build list of states we found
##    for trans in data:
##        transitions.append(trans['status_to'])
##
##    # build list of states we expect
##    should_have = ["ready","queued","pending","running","success"]
##
##    # test that we have all expected states
##    for state in should_have:
##        count = transitions.count(state)
##        if (count != 1):
##            sys.stderr.write("found unexpected number of states: found %s \"%s\" states" % (count,state))
##    #        assert 0
##        else:
##            print(f"found: {count} {state}")
##
##    # test that we have all expected states
##    cmd = 'jaws run task-status ' + run_id
##    rc,stdout,stderr = run(cmd)
##    data = json.loads(stdout)
##    
#    # testcase-2d
##    A list of tasks with unique IDs for the run jaws run task-status
#
#    # testcase-2e
#    if (not output_dir):
#        print("No output directory displayed in status output")
#        assert 0
#
#
#def test_task_level_info(submit_wdl):
#    output_dir = submit_wdl('output_dir')
#
#    # testcase-3a
#    # has state info & when in said state `task-status`
#    cmd = 'jaws run task-status ' + run_id
#    rc,stdout,stderr = run(cmd)
#    data = json.loads(stdout)
#
#    # testcase-3b
#    # Information since when task is in said state jaws run task-status
#
#    # testcase-3c
#    # A list of transitions between states (entered state at time, left state at time)  jaws run task-log
#
#    # testcase-3d
#    # stderr/stdout of the task (in future release)
#
#    # testcase-3e
#    # Shell script as executed by JAWS (in future release)
#
#
##def test_bad_input_path(run):
#    # testcase-4
#    # User should get helpful error message if they enter Invalid Input file like 
#    # bad path to input json or file is not json 
#
##def test_bad_input_json(run):
#    # testcase-5
#    # User should get helpful error message if they enter non existing paths in input file
#
##def test_bad_input_perms(run):
#    # testcase-6
#    # User should get helpful error message if they enter Inaccessible paths in input file (bad perms)
#
##def test_bad_wdl_syntax(run):
#    # testcase-7
#    # User should get helpful error message if they enter Invalid WDL (syntax) that womtool.jar will catch
#
##def test_bad_wdl_semantics(run):
#    # testcase-8
#    # User should get helpful error message if they enter Invalid WDL (semantics) cromwell.log error
#    # use wrong variable for input to task 2 for referencing_db_and_shifter
#
#
#def test_output_dir(submit_wdl):
#    output_dir = submit_wdl['output_dir']
#
#    # testcase-12a
#    # Output files, as specified during submission
#    if os.path.exists(output_dir) and os.path.isdir(output_dir):
#        if not os.listdir(output_dir):
#            print("Output directory is empty. It should have at least logs dir")
#    else:
#        print("Output directory doesn't exist")
#
#    # testcase-12b
#    # stderr/stdout of each task in the run
#
##def test_runtime_shared(run):
#    # testcase-20 (2.1 or later)
#    # shared: Test that runtime {shared: 1} works. 
#    # The shared: 1 option allows jaws user to re-use pools when they run another workflow and pools are still up.  
#    # The shared: 0 option only allows tasks within a run instance to share the pool.
#    # One thing to note is that if the workers removed by scancel, jtm doesn’t know it’s removed.
#    # once remove-pool is added to jaws command, you can use the command to remove workers.
#    # Test by submitting two of the same wdls with shared: 1 and verify only 1 pool gets created.  
#    # Then verify two pools get created with shared: 0.
#
#
##def test_runtime_constraint(run):
#    # testcase-21 (2.1 or later)
#    #constraint: what happens if you set constraint: "knl" and run on lbnl?  NERSC options are knl and haswell. LBNL has only haswell.
#    # What happens if you use an invalid constraint on cori like "jojo"? (hint: value get's ignored).
#
#
##def test_runtime_mem(run):
#    # TESTCASE-22 (2.1 or later)
#    # mem: what if this is out of bounds for the possible site resources. NERSC has limit of 450G machines.
#    # Test shared: 1 and mem: "600G" & there was already a pool with same name running. => job completed w/o error, which is not the desired outcome.
#    # Test shared: 0 and mem: "600G" & no pool running. => job hangs in running state and no job gets put into queue.
#
##def test_runtime_time(run):
#    # time: what if this is out of bounds for the possible site resources. NERSC has limit of 72hrs.
#    # I tested time: "99:99:99" (no pools running) and got error "Failed to get a reply from the manager:" 
#    # which is good; however, a better error message (and job should quite asap i.e. at client level).  
#    # Job was NOT stuck in "hanging" but said "failed".  (Warning: This error may have been due to 
#    # something else and not specific to this change).
#
