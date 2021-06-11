#!/usr/bin/env python
"""This script contains a function to submit multiple JAWS jobs,
and another that waits for all runs in a list to finish."""

import subprocess
import json
import time

def verify(runs, num_tasks):
    # check the number of rows in the task-log output for each task
    PASSED_RUN_IDS = []
    FAILED_RUN_IDS = []
    for run in runs:
        data = subprocess.run(['jaws', 'run', 'task-log', run], capture_output=True, text=True)
        taskLogOutput = data.stdout
        lines = taskLogOutput.splitlines()
        print(run + " has " + str(len(lines)) + " lines")

        # each task should have 5 transitions recorded and there is a header row
        if len(lines) == ((num_tasks * 5) + 1):
            print(run + " passes")
            PASSED_RUN_IDS.append(run)
        else:
            print(run + " fails")
            FAILED_RUN_IDS.append(run)

    print("Passed run ids: " + str(PASSED_RUN_IDS))
    print("Failed run ids: " + str(FAILED_RUN_IDS))
    print("passed: " + str(len(PASSED_RUN_IDS)))
    print("failed: " + str(len(FAILED_RUN_IDS)))
