#!/usr/bin/env python
"""This script reads a json inputs file to submit one or more wdls for automated testing.
   It appends a timestamp to the end of the output directory so that the same inputs file can be reused without
   overwriting existing data.
   It returns a dictionary version of the input file that has the run-id added for each input and the output directory
   updated with the timestamp info."""

import argparse
import sys
import os
import json
import parsing_functions as pf

def main():
    # parse arguments
    parser = argparse.ArgumentParser(
        description='This will submit one or more WDLs to JAWS and wait for completions.')
    parser.add_argument("-e", "--environment", help="The environement you want the job submitted to dev|staging|prod",
                        type=str, choices=['dev', 'staging', 'prod'], required=True)
    parser.add_argument("-i", "--inputs", help="The input json file that a list of wdls, input files and sites",
                        type=str, required=True)
    args = parser.parse_args()

    if not os.path.exists(args.inputs):
        print (r"Error: {args.inputs} does not exist")
        sys.exit(1)

    # source the environment specified on command line
    pf.source_environment(args.environment)

    with open(args.inputs) as json_file:
        # read the inputs file
        data = json.load(json_file)

        # add the environment specified on command line into the json data
        data['env'] = args.environment

        # submit each wdl
        for job in data['input-wdls']:
            # print('Site: ' + job['site'])

            run_id, time_dir = pf.submit_one_run(job['wdl'], job['inputs'], job['output-dir'], job['site'])

            # add the run_id and the updated output directory string to the dictionary info
            job['run-id'] = run_id
            job['output-dir'] = time_dir
            #run_id += 1  #delete this later

    print(json.dumps(data, indent=4))

if __name__ == '__main__':
    main()
