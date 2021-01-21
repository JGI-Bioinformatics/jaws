#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completion, 
   run tests and finally create a yaml file for the autoqc GUI."""

import sys
import os
import argparse
import json
import random
import string
import logging
import parsing_functions as pf

logging.basicConfig(filename='test_jaws_cmds.log', filemode='w', format='**%(asctime)s**\n%(message)s', level=logging.DEBUG)

# this is the name of the analysis file that will be created in the run's output directory
ANALYSIS_FILE_NAME = 'analysis.yaml'

# this must match the name of a threshold file that has been submitted to autoqc before running this test
THRESHOLD_FILE_NAME = 'test_jaws_cmds'

tmp_wdl = "pow23.wdl"
tmp_readme = "pow23.md"
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


def jaws_run_queue(final_dict,run_id,env):
    """ tests that the jaws run queue command has the run id in the stdout."""

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run queue | grep '\"id\":' | awk '{print $2}' | tr -d ','" % (env)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    ids=o.split()
    if str(run_id) in o:
        final_dict['queue'] = 1
    else:
        final_dict['queue'] = 0


def jaws_run_history(final_dict,run_id,env):
    """ tests that the jaws run history command has the run id in the stdout."""

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run history | grep '\"id\": %s' | awk '{print $2}' | tr -d ','" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    # if there was something in stdout, then grep found "id": <run_id>
    if o:
        final_dict['hist'] = 1
    else:
        final_dict['hist'] = 0


def jaws_wdl_add(final_dict,env):
    """ tests that the jaws wdl add command added something."""

    wdl= """workflow fq_count { 
         File fastq_file
         call count_seqs { input: infile = fastq_file }
         output { File outfile = count_seqs.outfile }   
        }

        task count_seqs {
        File infile
        command <<< echo ~{infile} >>> 
        output { File outfile = stdout() }
        }
        """

    # write a temporary wdl
    with open(tmp_wdl,"w") as f:
        f.write(wdl)

    # write a temporary readme
    readme='this workflow does not do anything'
    with open(tmp_readme,"w") as f:
        f.write(readme)

    # make sure this wdl doesn't already exist in the catalog (i.e. delete failed for the last test).
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list | grep %s" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    # if command succeeds, then there is still an old wdl in catalog.  We need to delete it.
    if r == 0:
        # delete wdl from catalog so we can test adding it back again
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl delete %s 2.0.0" % (env,wdl_catalog_name)
        (o,e,r) = pf.submit_cmd(cmd)
        logging.info("%s\n%s\n%s",cmd,o,e)
        if r:
            final_dict['add'] = "failed to delete old version"
            return
            
    # add to catalog
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl add %s 2.0.0 %s %s" % (env,wdl_catalog_name,tmp_wdl,tmp_readme)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    # show that it was added
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list | grep %s" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    # if there was something in stdout, then grep found "id": <run_id>
    if o:
        final_dict['add'] = 1
    else:
        final_dict['add'] = 0

def jaws_wdl_update(final_dict,env):
    """Update a readme for a WDL in the JAWS catalog"""
    readme = 'this readme has been changed'
    with open(tmp_readme,"w") as f:
        f.write(readme)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-doc %s 2.0.0 %s" % (env,wdl_catalog_name,tmp_readme)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)
    
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl about %s 2.0.0 " % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    if readme in o:
        final_dict['update'] = 1
    else:
        final_dict['update'] = 0

def jaws_wdl_update_wdl(final_dict,env):
    """Update a WDL from the JAWS catalog"""

    task = """task secondecho {
    command { echo second task }
    }"""

    with open(tmp_wdl,"a") as f:
        f.write(task)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl update-wdl %s 2.0.0 %s" % (env,wdl_catalog_name,tmp_wdl)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl get %s 2.0.0" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)
     
    if "secondecho" in o:
        final_dict['updateWDL'] = 1
    else:
        final_dict['updateWDL'] = 0

def jaws_wdl_versions(final_dict,env):
    """Test that we can the version for a given WDL
    "fq_count:1.0.0": {
        "created": "2020-11-04T22:24:11Z",
        "last_updated": "2020-11-04T22:24:11Z",
        "name": "fq_count",
        "owner": "jfroula",
        "production_release": "no",
        "version": "1.0.0"
        }
    """
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl versions %s | grep version | awk '{print $2}' " % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    if o.strip() == '\"2.0.0\"': 
        final_dict['versions'] = 1
    else:
        final_dict['versions'] = 0

def jaws_wdl_delete(final_dict,env):
    """Check that a WDL is deleted"""
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl delete %s 2.0.0" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws wdl list | grep %s" % (env,wdl_catalog_name)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    # if grep found wdl_catalog_name still in the catalog, delete failed
    if o:
        final_dict['delete'] = 0
    else:
        final_dict['delete'] = 1

def jaws_wdl_metadata(final_dict,run_id,env):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run metadata %s | grep workflowRoot | awk '{print $2}' | tr -d '\"'" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    if (o):
        final_dict['metadata'] = 1
    else:
        final_dict['metadata'] = 0

def jaws_wdl_errors(final_dict,run_id,env):
    """Check that a jaws run metadata returns workflowRoot has a value"""
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run errors %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    # we don't have any errors in our wdl so the errors command should't return anything.
    # We don't have a good test for this command here, but we only test that the return code is 0.
    if not r:
        final_dict['errors'] = 1
    else:
        final_dict['errors'] = 0
        
def jaws_wdl_task_status(final_dict,run_id,env):
    """Check that jaws run task-status returns something like this:
     fq_count.count_seqs 1   25177   running success 2021-01-13 12:37:45     The job completed successfully
    """
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-status %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    if 'fq_count.count_seqs' in o and 'The job completed successfully' in o:
        final_dict['task_status'] = 1
    else:
        final_dict['task_status'] = 0

def jaws_wdl_log(final_dict,run_id,env):
    """Check that the final line of jaws run log returns something like this:
        downloading download complete   2021-01-13 12:41:28 
    """
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run log %s | tail -1" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    if 'download complete' in o: 
        final_dict['log'] = 1
    else:
        final_dict['log'] = 0

def jaws_wdl_task_log(final_dict,run_id,env):
    """Check that the final line of jaws run log returns something like this:
        downloading download complete   2021-01-13 12:41:28 
    """
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run task-log %s" % (env,run_id)
    (o,e,r) = pf.submit_cmd(cmd)
    logging.info("%s\n%s\n%s",cmd,o,e)

    a=[]
    for line in o.split("\n"):
        a.append(line.split())

    # remove empty elements
    a = list(filter(None,a))

    try:
        if a[1][3] == 'created' and a[-1][4] == 'success': 
            final_dict['task_log'] = 1
    except:
        final_dict['task_log'] = 0

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

# uncomment for testing
#run_info = {'output_dir': '/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/1611111560207628', 'output_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'run_id': 16240, 'site_id': 'CORI', 'status': 'uploading', 'submission_id': 'e599f6b6-ed6a-4db5-8f0c-026a766a772c', 'upload_task_id': '8220363c-5acb-11eb-a4e5-0a53a3613b81'}

logging.info(f"submitted job: {run_info}\n")
run_id = run_info['run_id']

jaws_info(final_dict,args.environment)
jaws_status(final_dict,args.environment)
jaws_run_queue(final_dict,run_id,args.environment)

if pf.wait_for_one_run(args.environment,run_id,check_tries=check_tries,check_sleep=check_sleep):
    # if this function returns 1, that means 'jaws run status' worked. So lets set final_dict for status.
    final_dict['run_status'] = 1
else:
    sys.stderr.write("Workflow not complete after alloted time...exiting.")

jaws_run_history(final_dict,run_id,args.environment)
jaws_wdl_add(final_dict,args.environment)
jaws_wdl_update(final_dict,args.environment)
jaws_wdl_update_wdl(final_dict,args.environment)
jaws_wdl_versions(final_dict,args.environment)
jaws_wdl_delete(final_dict,args.environment)
jaws_wdl_metadata(final_dict,run_id,args.environment)
jaws_wdl_errors(final_dict,run_id,args.environment)
jaws_wdl_task_status(final_dict,run_id,args.environment)
jaws_wdl_log(final_dict,run_id,args.environment)
jaws_wdl_task_log(final_dict,run_id,args.environment)

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
logging.info("submit_analysis_file\n%s\n%s\n",o,e)
