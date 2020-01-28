"""
JAWS Analysis/Run management functions; these interact via REST with the JAWS Central server.
"""

import sys
import os
import pwd
import grp
import shutil
import json
import getpass
import requests
#import base64
#import html2text
import csv
import time
import shutil
import subprocess
from zipfile import ZipFile
import tarfile
import click
import uuid
import globus_sdk
from jaws_client import user, catalog, workflow

DEBUG = False
if "JAWS_DEBUG" in os.environ: DEBUG = True

# JAWS CONFIG
JAWS_URL = os.environ["JAWS_URL"]
jaws_workflows = JAWS_URL + "/wdl"
jaws_runs = JAWS_URL + "/run"
GLOBUS_ENDPOINT = os.environ["GLOBUS_ENDPOINT"]
GLOBUS_BASEDIR = os.environ["GLOBUS_BASEDIR"]
JAWS_STAGING_SUBDIR = os.environ["JAWS_STAGING_SUBDIR"]
JAWS_SITE = os.environ["JAWS_SITE"]

@click.group(context_settings={'help_option_names':['-h','--help']})
def run():
    """
    Run management.
    """
    pass

@run.command()
def queue():
    """
    List your unfinished runs.
    """
    current_user = user.User()
    try:
        r = requests.get(jaws_runs, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@run.command()
@click.option('--days', default=1)
def history(days, name=None):
    """
    List your past runs.
    """
    data = {
        "delta_days" : days
    }
    url = JAWS_URL + "/search"
    current_user = user.User()
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))
#    runs = []
#    if name:
#        for a_run in result['results']:
#            if 'name' in a_run and a_run['name'] == name:
#                runs.append(a_run)
#    else:
#        runs = result['results']
#    print(json.dumps(runs, indent=4, sort_keys=True))

def _run_status(run_id):
    """
    Return the status of a run in JSON format.
    """
    url = "%s/%s" % (jaws_runs, run_id)
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    return result

@run.command()
@click.argument('run_id')
def status(run_id):
    """
    Show current status of a run.
    """
    result = _run_status(run_id)
    print(json.dumps(result, indent=4, sort_keys=True))

@run.command()
@click.argument('run_id')
def tasks(run_id):
    """
    Show status of each task of a run.
    """
    metadata = _run_metadata(run_id)
    status = metadata["status"]
    print("status: " + status)

    # failures
    if "failures" in metadata:
        print("failures:")
        for failure in metadata["failures"]:
            print("\t" + failure["message"])

    # call summary
    print("calls:")
    calls = metadata['calls']
    result = []
    for task_name in calls.keys():
        task = calls[task_name]
        for attempt in task:
            status = attempt['executionStatus']
            shard = attempt['shardIndex']
            result.append((task_name, shard, status))
            if shard > -1:
                print(task_name + "[" + str(shard) + "]\t" + status)
            else:
                print("\t".join((task_name, status)))
            if status == "Failed":
                for failure in attempt["failures"]:
                    print("\t" + failure["message"])

@run.command()
@click.argument('run_id')
def wait(run_id):
    """
    Wait until run is complete; check return code.
    """
    while True:
        result = _run_status(run_id)
        if result['status'] == 'Succeeded':
            sys.exit(0)
        elif result['status'] == 'Failed' or result['status'] == 'Aborted':
            sys.exit(1)
        time.sleep(60)

def _run_metadata(run_id):
    """
    Return the metadata of a run.
    """
    url = "%s/%s/metadata" % (jaws_runs, run_id)
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    return result

@run.command()
@click.argument('run_id')
def metadata(run_id):
    """
    Return the detailed metadata for a run.
    """
    result = _run_metadata(run_id)
    print(json.dumps(result, indent=4, sort_keys=True))

#@run.command()
#@click.argument('run_id')
#def log(run_id):
#    """
#    Return the log of a run.
#    """
#    current_user = user.User()
#    try:
#        r = requests.get(jaws_runs + '/' + run_id + '/log', headers=current_user.header())
#    except:
#        sys.exit("Unable to communicate with JAWS server")
#    if r.status_code != 200:
#        sys.exit(r.text)
#    print(r.text)

@run.command()
@click.argument('run_id')
def cancel(run_id):
    """
    Cancel a run.
    """
    url = "%s/%s" % (jaws_runs, run_id)
    current_user = user.User()
    try:
        r = requests.delete(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code == 200:
        print('run ' + run_id + " was canceled")
    else:
        sys.exit(r.text)

#@run.command()
#@click.option('--days', default=0)
#def search(days=0, status=[]):
#    """
#    Search for a run by label.
#    """
#    global jaws_runs
#    data = {}
#    if days > 0:
#        data['delta_days'] = days
#    if status is not None and len(status) > 0:
#        data['status'] = status
##        for label in labels:
##                (key,value)=label.split(":")
##                data[key]=value
#    current_user = user.User()
#    try:
#        r = requests.post(jaws_runs + "/query", data=data, headers=current_user.header())
#    except:
#        sys.exit("Unable to communicate with JAWS server")
#    if r.status_code != 200:
#        sys.exit(r.text)
#    result = r.json()
#    print(json.dumps(result, indent=4, sort_keys=True))

@run.command()
@click.argument('run_id')
@click.option('--task', default=None)
def delete(run_id, task):
    """
    Delete the output of a run or task to avoid caching.
    """
    url = None
    if task is not None:
        url = "%s/%s" % (jaws_runs, run_id)
    else:
        url = "%s/%s/%s" % (jaws_runs, run_id, task_name)
    current_user = user.User()
    try:
        r = requests.delete(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    result = r.json()
    if r.status_code != 200:
        sys.exit(r.text)
    print(json.dumps(result, indent=4, sort_keys=True))

def _validate_workflow_name(workflow_full_name):
    """
    Set version to "latest" if not defined
    """
    if ":" not in args.workflow_full_name:
        workflow_full_name = workflow_full_name + ":latest"
    return workflow_full_name

def _labels_string_to_dict(labels_string):
    """
    Converts string to dictionary
    """
    labels_dict = {}
    labels = labels_string.split(',')
    for label in labels:
        if label.count(':') != 1:
            sys.exit("Invalid label: " + label + "; labels must be in \"key:value\" format")
        (key, value)=label.split(":")
        labels_dict[key] = value
    return labels_dict

#@run.command()
#@click.argument('submission_id')
#def block(submission_id):
#    """
#    Block until run is complete.
#    """
#    success = {
#        'Succeeded' : None
#    }
#    failure = {
#        'Failed' : None,
#        'Aborted' : None
#    }
#    while True:
#        result = _run_status(submission_id)
#        if result['status'] in success:
#            sys.exit(0)
#        elif result['status'] in failure:
#            sys.exit(1)
#        time.sleep(30)

         
# TODO add option for labels (JSON file?  kvargs?)
@run.command()
@click.argument('wdl_file', nargs=1)
@click.argument('infile', nargs=1)
@click.option('--outdir', default=None, help="Output dir")
@click.option('--site', default=None, help="Computing site")
#@click.option('--block', default=False, help="Block until complete")
def submit(wdl_file, infile, outdir, site): #, block):
    """
    Submit a run for execution.
    """
    current_user = user.User()

    # DEFINE OUTPUT DIR AND VERIFY IT'S ACCESSIBLE VIA GLOBUS
    if not outdir:
        if current_user.config["PATHS"]["default_outdir"]:
            outdir = current_user.config["PATHS"]["default_outdir"]
        else:
            sys.exit("--outdir required as no default specified in your config file\nYou may set your \"default_outdir\" by editing %s" % (current_user.config_file,))
    if not GLOBUS_BASEDIR: sys.exit("ERROR: \$GLOBUS_BASEDIR env var not defined")
    if not outdir.startswith(GLOBUS_BASEDIR): sys.exit("Invalid outdir, %s; Globus can only write under %s" % (outdir, GLOBUS_BASEDIR))

    # CREATE UNIQUE STAGING ID
    submission_uuid = str(uuid.uuid4())

    # VALIDATE RUN
    run = workflow.Workflow(wdl_file, infile)
    if not run.validate(): sys.exit("Failed validation; aborting")

    # PREPARE RUN
    staging_dir = os.path.join(GLOBUS_BASEDIR, JAWS_STAGING_SUBDIR)
    if not os.path.isdir(staging_dir): os.makedirs(staging_dir)
    run.prepare_wdls(staging_dir, submission_uuid)
    run.prepare_inputs(GLOBUS_BASEDIR, JAWS_STAGING_SUBDIR, JAWS_SITE, submission_uuid)

    # GET COMPUTE SITE INFO (E.G. GLOBUS ENDPOINT PARAMS)
    if not site: site = pick_site(run)
    dest = site_info(site)
    dest_dir = os.path.join(dest["basepath"], dest["staging_subdir"])

    # GLOBUS TRANSFER
    transfer_client = current_user.transfer_client()
    tdata = globus_sdk.TransferData(
        transfer_client,
        GLOBUS_ENDPOINT,
        dest["endpoint"],
        label=submission_uuid,
        sync_level="checksum",
        verify_checksum=True,
        preserve_timestamp=True,
        notify_on_succeeded=False,
        notify_on_failed=True,
        notify_on_inactive=True,
        skip_activation_check=False )
    for source_file, dest_file in run.manifest:
        abs_dest_file = os.path.join(dest_dir, dest_file)
        # NOTE: recursive=False means folders will not be transferred;
        # TODO: change this is WDL supports folders in the future.
        tdata.add_item(source_file, abs_dest_file, recursive=False)
    transfer_result = transfer_client.submit_transfer(tdata)
    transfer_task_id = transfer_result["task_id"]
    print("Staging files; Globus task_id = %s" % (transfer_task_id,))

    # SUBMIT RUN
    # NOTE THAT THE FILE TRANSFER IS NOT COMPLETE YET
    data = {
        "submission_uuid" : submission_uuid,
        "globus_endpoint" : GLOBUS_ENDPOINT,
        "outdir" : outdir,
        "globus_transfer_task_id" : transfer_task_id
    }
    try:
        r = requests.post(jaws_runs, data=data, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != requests.codes.ok:
        sys.exit(r.text)
    result = r.json()
    run_id = result["run_id"]
    print("Successfully queued run %s" % (run_id,))

    # OPTIONALLY WAIT UNTIL RUN IS COMPLETE
    #if block is True: block(submission_id)


def pick_site(run):
    """
    Ask Central to pick a computing Site based upon run attributes (e.g. RAM required, size of input data) and sites' availability.
    """
    print("Resources: RAM = %s MB; XFER = %s GB" % (run.max_ram_mb, run.transfer_gb))
    data = {
        "max_ram_mb" : run.max_ram_mb,
        "transfer_gb" : run.transfer_gb
    }
    url = "%s/run/pick_site" % (JAWS_URL,)
    current_user = user.User()
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    result = r.json()
    if r.status_code != 200: sys.exit(r.text)
    site = result["site"]
    return site


def get_site_info(site):
    """
    Get computing Site Globus endpoint parameters.
    """
    url = "%s/run/site_info" % (JAWS_URL,)
    current_user = user.User()
    try:
        r = requests.post(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    result = r.json()
    if r.status_code != 200: sys.exit(r.text)
    return result


@run.command()
def convert(batch_input_file, outfile=None):
    """
    Covert table of inputs to JSON format for submitting a batch of runs.
    """
    # read inputs table
    data = open(batch_input_file, "r")
    runs = csv.DictReader(data, delimiter="\t", quotechar="\"")
    batch = []
    line_num = 0
    for row in runs:
        # populate dictionary
        inputs = {}
        line_num = line_num + 1
        for key in runs.fieldnames:
            # validate table
            if not row[key]:
                sys.exit("Input table missing value for " + key + " in row " + str(line_num))
            inputs[key] = row[key]
        batch.append(inputs)

    # output
    if outfile:
        with open(batch_output_file, 'w') as outfile:
                json.dump(batch, outfile, indent=4, sort_keys=True)
        outfile.close()
    else:
        print(json.dumps(batch, indent=4, sort_keys=True))
