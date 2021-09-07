#!/usr/bin/env python
"""This script contains a function to submit multiple JAWS jobs,
and another that waits for all runs in a list to finish."""

import subprocess
import json
import time
import functools

# flush the printstream so that the output can be redirected to a file to prevent 
# loss of test results if the ssh connection is lost
print = functools.partial(print, flush=True)

SUBMIT_SLEEP = 1
CHECK_SLEEP = 360
CHECK_TRIES = 100

def submit(num_submissions, wdl, inputs, site):
    """Submit jobs to jaws and return a list of run ids."""

    # list to hold the run ids that are created
    run_ids = []

    # loop to do the submissions
    i = 1
    while i <= num_submissions:
        data = subprocess.run(['jaws', 'submit', '--no-cache', wdl, inputs, site],
                              capture_output=True, text=True)
        output = data.stdout

        # fake output used for testing without actually submitting jobs
        #output = '{"run_id": 12, "output_dir": "/global/u1/a/akollmer/mywdls/simple1/out4"}'

        out_json = json.loads(output)

        run_id = str(out_json["run_id"])
        run_ids.append(run_id)
        print(f'Loop: {i}   run_id: {run_id}   site: {site}')

        i += 1
        time.sleep(SUBMIT_SLEEP)

    print(run_ids)
    assert len(run_ids) == num_submissions
    return run_ids


def wait(run_ids):
    SUCCESSFUL_RUN_IDS = []
    FAILED_RUN_IDS = []

    """Wait for all the runs in run_ids list to finish."""
    for run in run_ids:
        tries = 1
        while tries < CHECK_TRIES:
            # check whether the run has finished every 60 seconds
            data = subprocess.run(['jaws', 'status', run], capture_output=True, text=True)
            output = data.stdout

            # fake output used for testing without actually submitting jobs
            #output = '{"status": "download complete"}'

            status_json = json.loads(output)

            run_status = status_json["status"]

            print("Run " + run + " status after " + str(tries) + " try is: " + run_status)

            if run_status == "download complete":
                # check whether the run succeeded
                # "result": "succeeded",
                run_result = status_json["result"]
                if run_result == "succeeded":
                    print(run + " succeeded")
                    SUCCESSFUL_RUN_IDS.append(run)
                else:
                    print(run + " failed")
                    FAILED_RUN_IDS.append(run)
                break
            else:
                time.sleep(CHECK_SLEEP)
            tries += 1

    print("All runs have finished")
    print("Successful run ids: " + str(SUCCESSFUL_RUN_IDS))
    print("Failed run ids: " + str(FAILED_RUN_IDS))
    print("Successful count: " + str(len(SUCCESSFUL_RUN_IDS)))
    print("Failed count: " + str(len(FAILED_RUN_IDS)))
    print("===========================================")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        submit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        NUM_SUBMISSIONS = 6
        WDL = '/global/u1/a/akollmer/mywdls/simple3/simple3.wdl'
        INPUTS = '/global/u1/a/akollmer/mywdls/simple3/simple3_inputs.json'
        SITE = 'cori'
        submit(NUM_SUBMISSIONS, WDL, INPUTS, SITE)
