#!/usr/bin/env python
"""This script submits multiple JAWS jobs, waits for them to finish and then
checks that each run contains the correct number of entries in the task-log.
User must supply correct value for NUM_TASKS_IN_WDL for the wdl being used."""

import subprocess
import submit_loop

NUM_SUBMISSIONS = 20
WDL = '/global/u1/a/akollmer/mywdls/simple3/simple3.wdl'
NUM_TASKS_IN_WDL = 3
INPUTS = '/global/u1/a/akollmer/mywdls/simple3/simple3_inputs.json'
OUT_DIR_ROOT = '/global/cscratch1/sd/akollmer/loop/'
SITE = 'nersc'

# submit jobs to jaws
RUNS = submit_loop.submit(NUM_SUBMISSIONS, WDL, INPUTS, OUT_DIR_ROOT, SITE)

#wait for jobs to finish
submit_loop.wait(RUNS)

# check the number of rows in the task-log output for each task
PASSED_RUN_IDS = []
FAILED_RUN_IDS = []
for run in RUNS:
    data = subprocess.run(['jaws', 'run', 'task-log', run], capture_output=True, text=True)
    taskLogOutput = data.stdout
    lines = taskLogOutput.splitlines()
    print(run + " has " + str(len(lines)) + " lines")

    # each task should have 5 transitions recorded and there is a header row
    if len(lines) == ((NUM_TASKS_IN_WDL * 5) + 1):
        print(run + " passes")
        PASSED_RUN_IDS.append(run)
    else:
        print(run + " fails")
        FAILED_RUN_IDS.append(run)

print("Passed run ids: " + str(PASSED_RUN_IDS))
print("Failed run ids: " + str(FAILED_RUN_IDS))
print("passed: " + str(len(PASSED_RUN_IDS)))
print("failed: " + str(len(FAILED_RUN_IDS)))
