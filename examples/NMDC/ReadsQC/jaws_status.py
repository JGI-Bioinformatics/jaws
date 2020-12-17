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


def parse_submit_log(filename):
    run_ids = []
    with open(filename, 'r') as fh:
        for line in fh:
            outs = []
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split()
            outs.append(cols[0])
            if len(cols) > 1:
                outs.append(cols[1])
            run_ids.append(outs)
    return run_ids

    
def main():
    if len(sys.argv) != 2:
        print(f"Usage: {__file__} submit.log")
        sys.exit(1)

    print("# RunId\tStatus\tResult\tCromwellId\tInput")
    for ref in parse_submit_log(sys.argv[1]):
        run_id = ref[0]
        label = ref[1] if len(ref) > 1 else ''
        json_data = run_jaws_status(run_id)
        print("%s\t%s\t%s\t%s\t%s"%(run_id, json_data.get('status', None), json_data.get('result', None), json_data.get('cromwell_run_id', None), label))


if __name__ ==  '__main__':
    main()
