#!/usr/bin/env python
'''
Helper script to submit jaws runs.
Inputs are wdl name, site nane (cori, jgi), file of files

The file of files is a list of input files for the inputs.json.
'''

import os
import sys
import subprocess
import tempfile
import json
import shutil

def run_jaws(wdl, json_file, site):
    cmd = f"jaws run submit {wdl} {json_file} {site}"
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return str(stdout.decode('utf-8')).strip(), str(stderr.decode('utf-8')).strip(), process.returncode


def get_run_id(stdout):
    json_str = ''
    for line in stdout.split('\n'):
        if line.startswith('{'):
            json_str += line.strip()
        elif json_str:
            if line.endswith('}'):
                json_str += line.strip()
                break
            else:
                json_str += line.strip()
    json_data = json.loads(json_str)
    return json_data.get("run_id", None)


def get_input_files(filename):
    files = []
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            files.append(line)
    return files


def create_input_json(json_data, json_file):
    with open(json_file, 'w') as fh:
        json.dump(json_data, fh)


def main():
    if len(sys.argv) != 4:
        print("Usage %s <wdl> <site> <file_of_input_files>"%__file__)
        sys.exit()

    valid_sites = ['cori', 'jgi']

    wdl = sys.argv[1]
    site = sys.argv[2]
    fof = sys.argv[3]
    run_ids = []
    submit_log_file = "submit.log"
    tmp_dir = tempfile.mkdtemp()

    if site.lower() not in valid_sites:
        raise SystemExit("Invalid site input.")

    ## Need to change this to match inputs.json for specified wdl.
    input_json_data = {
        "jgi_rqcfilter.input_files": [
        ],
        "jgi_rqcfilter.database": "/refdata/nmdc/RQCFilterData"
    }

    print(f"Creating {submit_log_file}")
    fh = open(submit_log_file, "w")

    infiles = get_input_files(fof)
    for infile in infiles:
        # Specific inputs.json entry for specified wdl.
        input_json_data["jgi_rqcfilter.input_files"] = [infile]

        # create tmp inputs.json file
        input_json_file = os.path.join(tmp_dir, "%s.json"%os.path.basename(infile))
        create_input_json(input_json_data, input_json_file)

        # Submit jaws run
        stdout, stderr, exitcode = run_jaws(wdl, input_json_file, site)
        run_id = get_run_id(stdout)
        run_ids.append([run_id, infile])

        fh.write(f"{run_id}\t{infile}\n")

    # remove temp inputs json files
    shutil.rmtree(tmp_dir)

    fh.close()


if __name__ == '__main__':
    main()
