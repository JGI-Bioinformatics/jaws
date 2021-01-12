#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completions.
   It will then submit a test cmd and create an analysis.yaml file."""

import sys
import os
import argparse
import json
import parsing_functions as pf

# this is the name of the analysis file that will be created in the run's output directory
ANALYSIS_FILE_NAME = 'analysis.yaml'
# this must match the name of a threshold file that has been submitted to autoqc before running this test
THRESHOLD_FILE_NAME = 'jaws_run_success'

# parse arguments
parser = argparse.ArgumentParser(
    description="This will submit a WDL to JAWS and wait for completions. It will then submit "
                "a test cmd and create an analysis.yaml file to be used for the auto-qc command.")
parser.add_argument("-w", "--wdl", help="The WDL you want to submit to JAWS", type=str, required=True)
parser.add_argument("-i", "--inputs", help="The input json file that accompanies the WDL", type=str, required=True)
parser.add_argument("-s", "--site", help="The site at which you wish to run [cori|jgi]", type=str, required=True)
parser.add_argument("-e", "--environment", help="The JAWS environment to run the test in [prod|staging|dev]",
                    type=str, required=True)
args = parser.parse_args()

if not os.path.exists(args.wdl):
    print(r"Error: {args.wdl} does not exist")
    sys.exit(1)

if args.environment.lower() not in ['prod', 'staging', 'dev']:
    print(r"Error: environment must be prod, staging or dev")

####################
#     MAIN         #
####################
#
# Submit the WDL to JAWS and
# wait for jaws run to complete
#
run_json = pf.submit_one_run_to_env(args.wdl, args.inputs, ".", args.site, args.environment)
run_info = json.loads(run_json)
run_id = run_info['run_id']

pf.wait_for_one_run(run_id)

#
# Run: jaws run status
#
cmd = "jaws run status %s" % run_id
(o, e, r) = pf.submit_jaws_cmd(cmd, args.environment)
status_info = json.loads(o)

# Write the status and result of the run into a dictionary
# that will become the AutoQC thresholds file as yaml.
final_dict = {}
final_dict['status'] = status_info['status']
final_dict['result'] = status_info['result']

#
# create the analysis yaml file
#
# put the analysis file in the run's output folder so that running multiple tests
# at one time does not result in the analysis files overwriting one another
analysis_file_path = status_info['output_dir'] + '/' + ANALYSIS_FILE_NAME

# create the name that will show up in the AutoQC report header by appending
# the test name 'run_success' to the wdl file name
test_report_name = os.path.basename(args.wdl) + '_run_success'
pf.create_analysis_file(final_dict, analysis_file_path, args.name)

#
# submit the analysis.yml file
#
pf.submit_analysis_file(analysis_file_path, THRESHOLD_FILE_NAME, args.environment)
