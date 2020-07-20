#!/usr/bin/env python
"""This script contains a function to submit multiple JAWS jobs,
and another that waits for all runs in a list to finish."""

import subprocess
import json
import time

SUBMIT_SLEEP = 1
CHECK_SLEEP = 60
CHECK_TRIES = 100

def submit(num_submissions, wdl, inputs, out_dir_root, site):
    """Submit jobs to jaws and return a list of run ids."""

    # create timestamp string to make output directory unique
    ts_str = str(time.time())
    ts_str = ts_str.replace(".", "")
    print(ts_str)

    # list to hold the run ids that are created
    run_ids = []

    # loop to do the submissions
    i = 1
    while i <= num_submissions:
        outdir = out_dir_root + ts_str + "/out" + str(i)
        data = subprocess.run(['jaws', 'run', 'submit', wdl, inputs, outdir, site],
                              capture_output=True, text=True)
        output = data.stdout

        # fake output used for testing without actually submitting jobs
        #output = '{"run_id": 12, "output_dir": "/global/u1/a/akollmer/mywdls/simple1/out4"}'

        out_json = json.loads(output)

        run_id = str(out_json["run_id"])
        out_dir = str(out_json["output_dir"])
        run_ids.append(run_id)
        print(f'Loop: {i}   run_id: {run_id}   output_dir: {out_dir}')

        i += 1
        time.sleep(SUBMIT_SLEEP)

    print(run_ids)
    assert len(run_ids) == num_submissions
    return run_ids


def wait(run_ids):
    """Wait for all the runs in run_ids list to finish."""
    for run in run_ids:
        tries = 1
        while tries < CHECK_TRIES:
            # check whether the run has finished every 60 seconds
            time.sleep(CHECK_SLEEP)
            data = subprocess.run(['jaws', 'run', 'status', run], capture_output=True, text=True)
            output = data.stdout

            # fake output used for testing without actually submitting jobs
            #output = '{"status": "download complete"}'

            status_json = json.loads(output)

            run_status = status_json["status"]
            print("Run " + run + " status after " + str(tries) + " try is: " + run_status)

            if run_status == "download complete":
                break

            tries += 1

    print("All runs have finished")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        submit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        NUM_SUBMISSIONS = 6
        WDL = '/global/u1/a/akollmer/mywdls/simple3/simple3.wdl'
        INPUTS = '/global/u1/a/akollmer/mywdls/simple3/simple3_inputs.json'
        OUT_DIR_ROOT = '/global/cscratch1/sd/akollmer/loop/'
        SITE = 'nersc'
        submit(NUM_SUBMISSIONS, WDL, INPUTS, OUT_DIR_ROOT, SITE)
