#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completion, 
   run tests and finally create a yaml file for the autoqc GUI."""

import sys
import os
import argparse
import json
import logging
import parsing_functions as pf

logging.basicConfig(filename='bad_inputs_msg.log',
                    filemode='w',
                    format='**%(asctime)s**\n%(message)s',
                    level=logging.DEBUG)

# this is the name of the analysis file that will be created in the run's output directory
ANALYSIS_FILE_NAME = 'analysis.yaml'

# this must match the name of a threshold file that has been submitted to autoqc before running this test
THRESHOLD_FILE_NAME = 'jaws_run_success'

check_tries = 50  # try this many times when waiting for a JAWS run to complete.
check_sleep = 30  # wait for this amount of time between tries.

#########################
#    Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#

def json_file_does_not_exist(e):
    # check wh
    if "File(s) not accessible:" in e:
        final_dict['json_missing'] = 1
    else:
        final_dict['json_missing'] = 0

def test_json_bad_path_to_input_file_msg():
    pass


####################
#     MAIN
####################

# parse arguments
parser = argparse.ArgumentParser(
description="This will submit a WDL to JAWS and wait for completions. It will then submit "
            "a test cmd and create an analysis.yaml file to be used for the auto-qc command.")
parser.add_argument("-s", "--site", help="The site at which you wish to run [cori|jgi]", type=str, required=True)
parser.add_argument("-e", "--environment", help="The JAWS environment to run the test in [prod|staging|dev]",
                    type=str, required=True)
args = parser.parse_args()

if args.environment.lower() not in ['prod', 'staging', 'dev']:
    print(r"Error: environment must be prod, staging or dev")

final_dict = {}

##########################
# Submit WDLs & run Tests #
##########################
#
# Submit path to json file that does not exist
# Can't use pf.submit_one_run_to_env here because it exits if submission not successful
source_cmd = "source ~/jaws-%s.sh > /dev/null && " % args.environment
wdl = "../TestsWDLs/jaws-alignment-example/main.wdl"
json = "./FileDoesNotExist.json"
out_dir = "out"  # nothing will get written here because these submissions are not accepted

submit_cmd = "jaws run submit %s %s %s %s" % (wdl, json, out_dir, args.site)

cmd = source_cmd + submit_cmd
(o, e, r) = pf.submit_cmd(cmd)
logging.info("cmd: %s\nout: %s\nerror: %s", cmd, o, e)

json_file_does_not_exist(e)


# Submit json that points to non-existent input file
# Can't use pf.submit_one_run_to_env here because it exits if submission not successful
# json = "../TestsWDLs/bad-inputs/bad_path.json"
# submit_cmd = "jaws run submit %s %s %s %s" % (wdl, json, out_dir, args.site)


#################################
# create the analysis yaml file #
#################################
#
# Put the analysis file in the run's output folder so that running multiple tests
# at one time does not result in the analysis files overwriting one another
#analysis_file_path = run_info['output_dir'] + '/' + ANALYSIS_FILE_NAME

# create the name that will show up in the AutoQC report header by appending
# the test name 'run_success' to the wdl file name
#test_report_name = os.path.basename(args.wdl) + '_' + str(run_id)
#logging.info(f"test_report_name: {test_report_name}\n")

#pf.create_analysis_file(final_dict, analysis_file_path, test_report_name)
#logging.info("create_analysis_file %s\n",analysis_file_path)

# submit the analysis.yml file
#(o, e, r) = pf.submit_analysis_file(analysis_file_path, THRESHOLD_FILE_NAME, args.environment)
#logging.info("submit_analysis_file\n%s\n%s", o, e)
