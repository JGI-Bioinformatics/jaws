#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completions.
   It will then submit a test cmd and create an analysis.yaml file."""

import subprocess
import sys
import os
import argparse
import tempfile
import json
import parsing_functions as pf

# parse arguments
parser = argparse.ArgumentParser(description='This will submit a WDL to JAWS and wait for completions.  It will then submit a test cmd and create an analysis.yaml file to be used for the auto-qc command.')
parser.add_argument("-w","--wdl",help="The WDL you want to submit to JAWS", type=str, required=True)
parser.add_argument("-i","--inputs",help="The input json file that accompanies the WDL", type=str, required=True)
parser.add_argument("-a","--analysis",help="The name you want to call the analysis yaml file", type=str, required=True)
parser.add_argument("-s","--site",help="The site at which you wish to run [cori|jgi]", type=str, required=True)
parser.add_argument("-n","--name",help="The name that will show up in the AutoQC report header", type=str, required=True)
args = parser.parse_args()

if not os.path.exists(args.wdl):
    print (r"Error: {args.wdl} does not exist")
    sys.exit(1)

####################
###     MAIN     ###
####################
#
# Submit the WDL to JAWS and
# wait for jaws run to complete
#
run_id = pf.submit_one_run(args.wdl,args.inputs,".",args.site)
pf.wait_for_one_run(run_id)

#
# Run: jaws run status
#
cmd = "jaws run status %s" % run_id
(o,e,r) = pf.submit_cmd(cmd)
data = json.loads(o)

# Varify output from the jaws command looks right.
#
# This dictionary will become the AutoQC thresholds file as yaml.
# It needs to be global to each test since each test will add a key/value pair.
final_dict = {}
final_dict['status'] = data['status']
final_dict['result'] = data['result']

#
# create the analysis yaml file
#
pf.create_analysis_file(final_dict,args.analysis,args.name)
