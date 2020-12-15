#!/usr/bin/env python
"""
Helper script to report the status of jaws runs by run_id.

Usage: jaws_status.py <start_run_id> <end_run_id>

Prints the run id, status, result and cromwell run id for the range
of input run ids
"""

import os
import sys
import subprocess
import json

def run_jaws_status(run_id):
    json_data = {}
    cmd = f"jaws run status {run_id}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if not process.returncode and stdout:
        json_data = json.loads(stdout)
    return json_data


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"Usage: {__file__} <start_run_id> <end_run_id>")

    run_id_start = int(sys.argv[1])
    run_id_stop = int(sys.argv[2])
    print("# RunId\tStatus\tResult\tCromwellId")
    for run_id in range(run_id_start, run_id_stop+1):
        json_data = run_jaws_status(run_id)
        print("%s\t%s\t%s\t%s"%(run_id, json_data['status'], json_data['result'], json_data['cromwell_run_id']))


if __name__ ==  '__main__':
    main()
