#!/usr/bin/env python
"""This script contains functions to submit jaws commands to test and functions to parse the output of these submissions"""

import subprocess
import json
import time
import sys
import logging

SUBMIT_SLEEP = 2

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


def submit_one_run_to_env(wdl, inputs, dir, site, env):
    """Submit a single job to jaws and return the run info as a dictionary."""

    # create timestamp string to make output directory unique
    out_dir = timestamp_dir(dir)

    # source the jaws environment and submit the run
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run submit %s %s %s %s" % (env, wdl, inputs, out_dir, site)
    (out, err, rc) = submit_cmd(cmd)

    if (rc > 0):
        print(f"Failed cmd: {cmd}\n{out}\n{err}\n")
        sys.exit(1)

    return json.loads(out)


def wait_for_multi_runs(run_ids,check_tries=100,check_sleep=30):
    """Wait for all the runs in run_ids list to finish."""
    for run in run_ids:
        tries = 1
        while tries < check_tries:
            # check whether the run has finished every 60 seconds
            time.sleep(check_sleep)
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

def wait_for_one_run(env,run_id,check_tries=100,check_sleep=30):
    """Wait for all the runs in run_ids list to finish."""
    run_id = str(run_id)
    tries = 1
    while tries <= check_tries:
        # check whether the run has finished every 60 seconds
        time.sleep(check_sleep)
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env,run_id)
        (o,e,r) = submit_cmd(cmd)

        status_json = json.loads(o)

        run_status = status_json["status"]
        logging.info("Run " + run_id + " status after " + str(tries) + " try is: " + run_status)

        if run_status == "download complete":
            return 1

        tries += 1

    logging.info("All tries have have been exhausted")
    return 0

def submit_cmd(cmd):
    """returns the stdout, stderr and rc"""
    process=subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False)
    output = process.stdout
    stderror = process.stderr
    thereturncode = process.returncode
    #print(cmd)

    return (output,stderror,thereturncode)

def submit_analysis_file(analysis_file, threshold_file_name, env):
    #  build a command like below
    #    global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/bin/autoqc_tool_gp.py \
    #    -p jaws qc analysis.yaml jaws_status

    autoqc_script = 'source /global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/config/rqc38.sh && \
    /global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/bin/autoqc_tool_gp.py'

    # make string to specify autoqc db for the JAWS environment we are using
    # commented out for now because we only have one jaws autoqc db currently
    # autoqcDb = 'jaws' + "_" + env
    autoqcDb = 'jaws'  # for now just use jaws

    cmd = autoqc_script + " -p " + autoqcDb + " qc " + analysis_file + " " + threshold_file_name
    logging.info(cmd)
    return submit_cmd(cmd)


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
                sys.stderr.write("The value for the following dictionary value cannot be determined. It should be an int, str, or bool.\n%s: %s\n" % (key,value))
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

