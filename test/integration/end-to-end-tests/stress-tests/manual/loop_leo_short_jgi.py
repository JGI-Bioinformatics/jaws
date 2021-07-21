#!/usr/bin/env python
"""This script submits multiple JAWS jobs, waits for them to finish and then
checks that each run contains the correct number of entries in the task-log.
User must supply correct value for NUM_TASKS_IN_WDL for the wdl being used."""

import subprocess
import submit_loop
import verify_task_log

NUM_SUBMISSIONS = 10
WDL = '/global/homes/a/akollmer/jaws/examples/leo_dapseq/Azospirillum_brasilense.wdl'
NUM_TASKS_IN_WDL = 15
INPUTS = '/global/homes/a/akollmer/jaws/examples/leo_dapseq/shortened.json'
SITE = 'jgi'
VERIFY_TASK_LOG_STATUS = False

# submit jobs to jaws
RUNS = submit_loop.submit(NUM_SUBMISSIONS, WDL, INPUTS, SITE)

#wait for jobs to finish
submit_loop.wait(RUNS)

# verify the task_log status
if VERIFY_TASK_LOG_STATUS:
    verify_task_log.verify(RUNS, NUM_TASKS_IN_WDL)
