#!/usr/bin/env python
"""This script contains functions to submit jaws commands to test and functions to parse the output of these submissions"""

import subprocess
import json
import time
import sys
import os
import pprint
import shlex

SUBMIT_SLEEP = 2
CHECK_SLEEP = 30
CHECK_TRIES = 100

def submit_multi_runs(num_submissions, wdl, inputs, dir, site):
    """Submit jobs to jaws and return a list of run ids."""

    # create timestamp string to make output directory unique
    out_dir_root = timestamp_dir(dir)
    print("output directory root: " + out_dir_root)

    # list to hold the run ids that are created
    run_ids = []

    # loop to do the submissions
    i = 1
    while i <= num_submissions:
        outdir = out_dir_root + "/out" + str(i)
        data = subprocess.run(['jaws', 'run', 'submit', wdl, inputs, outdir, site],
                              capture_output=True, text=True)
        output = data.stdout

        # fake output used for testing without actually submitting jobs
        #output = '{"run_id": 737, "output_dir": "/global/u1/a/akollmer/mywdls/simple1/out4"}'

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


def submit_one_run(wdl, inputs, dir, site):
    """Submit a job to jaws and return the run id and output directory."""

    # create timestamp string to make output directory unique
    out_dir = timestamp_dir(dir)
    print("output directory: " + out_dir)

    data = subprocess.run(['jaws', 'run', 'submit', wdl, inputs, out_dir, site],
                          capture_output=True, text=True)
    output = data.stdout

    # fake output used for testing without actually submitting jobs
    #output = '{"run_id": 741, "output_dir": "/global/cscratch1/sd/jfroula/JAWS/AutoQC/out"}'

    out_json = json.loads(output)

    run_id = str(out_json["run_id"])
    out_dir = str(out_json["output_dir"])
    print(f'run_id: {run_id} output_dir: {out_dir}')

    return run_id, out_dir


def wait_for_multi_runs(run_ids):
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

def wait_for_one_run(run_id):
    """Wait for all the runs in run_ids list to finish."""
    tries = 1
    while tries < CHECK_TRIES:
        # check whether the run has finished every 60 seconds
        time.sleep(CHECK_SLEEP)
        data = subprocess.run(['jaws', 'run', 'status', run_id], capture_output=True, text=True)
        output = data.stdout

        # fake output used for testing without actually submitting jobs
        #output = '{"status": "download complete"}'

        status_json = json.loads(output)

        run_status = status_json["status"]
        print("Run " + run_id + " status after " + str(tries) + " try is: " + run_status)

        if run_status == "download complete":
            break

        tries += 1

    print("All tries have have been exhausted")

def submit_cmd(cmd):
    process=subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False)
    output = process.stdout
    stderror = process.stderr
    thereturncode = process.returncode

    # combine stdout & stderr since they can sometimes be mixed up for some scripts.
    er=output + "\n" + stderror

    if thereturncode >= 1:
        sys.stderr.write("Error: command failed\n%s\n%s\n%s" % (cmd,stderror,output))
        sys.exit(1)
    else:
        print(f"command successfully submitted: {cmd}\n{er}")

    return (output,stderror,thereturncode)

def source_environment(env):
    cmd = 'source ~/jaws-{}.sh'.format(env)
    command = shlex.split("bash -c '{} && env'".format(cmd))
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    for line in proc.stdout:
        (key, _, value) = line.partition("=")
        os.environ[key] = value
    proc.communicate()

    pprint.pprint(dict(os.environ))
    return cmd

def create_analysis_file(final_dict,analysis_file,test_name):
    # Create the new analysis yaml file
    with open(analysis_file,"w") as wh:
        wh.write("metadata:\n")
        wh.write(f"    name: {test_name}\n")
        wh.write("data:\n")
        for key,value in final_dict.items():
            if (isinstance(value, (int, float, bool))):
                wh.write("    %s: %s\n" % (key,value))
            elif (isinstance(value, (str))):
                wh.write("    %s: \"%s\"\n" % (key,value))
            else:
                sys.stderr.write("The value of the input json file can not be determined. It should be an int, str, or bool.")
                sys.exit(1)

def timestamp_dir(dir):
    # create timestamp string to make output directory unique
    ts = str(time.time()).replace(".", "")
    if dir[-1] != '/':
        dir = dir + '/'
    return dir + ts

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        submit_multi_runs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

