#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completion, 
   run tests and finally create a yaml file for the autoqc GUI."""

import sys
import os
import argparse
import json
import logging
import parsing_functions as pf

logging.basicConfig(filename='test_jaws_cmds.log', filemode='w', format='**%(asctime)s**\n%(message)s', level=logging.DEBUG)

# this is the name of the analysis file that will be created in the run's output directory
ANALYSIS_FILE_NAME = 'analysis.yaml'

# this must match the name of a threshold file that has been submitted to autoqc before running this test
THRESHOLD_FILE_NAME = 'jaws_run_success'

wdl_catalog_name="tmp_wdl_catalog_name"
check_tries = 50 # try this many times when waiting for a JAWS run to complete.
check_sleep = 30 # wait for this amount of time between tries.

#########################
###     Functions     ###
#########################
#
# Test functions for verification of jaws log commands (log,task-log,status,task-status).
#
def jaws_info(final_dict,env):
    """ tests that there is a valid output for jaws info. Name should be dev,staging, or prod and version should have some value.
    {
    "docs_url": "https://jaws-docs.readthedocs.io/en/latest/",
        "name": "prod",
        "version": "2.1"
    }
    """

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws info" % (env)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    data = json.loads(o)

    # do we have an acceptable name
    if (data["name"] in ["prod","staging","dev"]) and data["version"] is not None:
        final_dict['info'] = 1
    else:
        final_dict['info'] = 0

def jaws_status(final_dict,env):
    """ tests that the jaws status is working. We don't care if some services are down.
        Just test that all below services are shown, regardless of status.
    {
        "CORI-Cromwell": "UP",
        "CORI-RMQ": "UP",
        "CORI-Site": "UP",
        "JAWS-Central": "UP",
        "JGI-Cromwell": "UP",
        "JGI-RMQ": "UP",
        "JGI-Site": "UP"
    }
    """

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status" % (env)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)
    data = json.loads(o)

    actual_keys = list(data.keys())
    required_keys = ["CORI-Cromwell","CORI-RMQ","CORI-Site","JAWS-Central","JGI-Cromwell","JGI-RMQ","JGI-Site"]

    bad=0
    for k in required_keys:
        if k not in actual_keys:
            bad=1

    if bad:
        # something bad happened
        final_dict['status'] = 0
    else:
        #success
        final_dict['status'] = 1


def jaws_wdl_metadata(final_dict,run_id,env):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run metadata %s | grep workflowRoot | awk '{print $2}' | tr -d '\"'" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    if (o):
        final_dict['metadata'] = 1
    else:
        final_dict['metadata'] = 0


####################
###     MAIN     ###
####################

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

final_dict = {}

##########################
# Submit WDL & run Tests #
##########################
#
# Submit the WDL to JAWS and
# wait for jaws run to complete
#
run_info = pf.submit_one_run_to_env(args.wdl, args.inputs, ".", args.site, args.environment)
run_id = run_info['run_id']

jaws_info(final_dict,args.environment)
jaws_status(final_dict,args.environment)

if not pf.wait_for_one_run(args.environment,run_id,check_tries=check_tries,check_sleep=check_sleep):
    sys.stderr.write("Workflow not complete after alloted time...exiting.")

jaws_wdl_metadata(final_dict,run_id,args.environment)

#print(final_dict)

#################################
# create the analysis yaml file #
#################################
#
# Put the analysis file in the run's output folder so that running multiple tests
# at one time does not result in the analysis files overwriting one another
analysis_file_path = run_info['output_dir'] + '/' + ANALYSIS_FILE_NAME

# create the name that will show up in the AutoQC report header by appending
# the test name 'run_success' to the wdl file name
test_report_name = os.path.basename(args.wdl) + '_' + str(run_id)
logging.info(f"test_report_name: {test_report_name}\n")

pf.create_analysis_file(final_dict, analysis_file_path, test_report_name)
logging.info("create_analysis_file %s\n",analysis_file_path)

# submit the analysis.yml file
(o,e,r) = pf.submit_analysis_file(analysis_file_path, THRESHOLD_FILE_NAME, args.environment)
logging.info("submit_analysis_file\n%s\n%s",o,e)
