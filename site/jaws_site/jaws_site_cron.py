#!/usr/bin/env python3

"""
JAWS Site needs to perform regularly scheduled tasks, such as post-processing of recently completed runs.  Each task is it's own function.  We are using this instead of cron, which isn't available on all systems and we want all systems to be maintained as nearly identically as possible.
"""

import sys
import schedule
import time
import os
import re
import subprocess
import pathlib
import requests
import json
#import pprint





#import config
#from jaws_site import config
#config = config.Config(env="JAWS_SITE_CONFIG")


#CROMWELL_URL = os.environ["CROMWELL_URL"]
GLOBUS_ENDPOINT =  # EP of this Site


def check_uploads():
    """
    Check Globus transfer status for uploads
    """
    # QUERY DB FOR ALL RUN WHICH ARE UPLOADING


def _check_upload():
    """
    Check on an uploading run.
    """
    # CHECK STATUS
    owner = user.User(user_id)
    transfer_client = owner.transfer_client()
    task = transfer_client.get_task(transfer_task_id)
    status = task["status"]
    if status == "ACTIVE": return
    elif status == "SUCCEEDED":
        # RUN IS READY, POST TO CROMWELL AND RECORD NEW STATE
        # TODO

    elif state == "FAILED":
        # UPDATE STATE

    #elif state == "INACTIVE": return

#############

def check_runs():
    """
    Check on runs that have been submitted to Cromwell.
    """
    # GET ALL RUNS WHICH HAVE BEEN SUBMITTED TO CROMWELL
    runs = [] # TODO
    for run in runs:
        _check_running(run)

def _check_run(user_id, run_id, cromwell_id, dest_endpoint, dest_dir):
    """
    Query Cromwell for current state
    """
    # GET STATUS FROM CROMWELL
    r = requests.get("%s/%s/status" % (self.cromwell, "run_id"))
    if r.status_code != requests.codes.ok: return
    status = r.json()["status"]

    if status == "Failed":
        # update db record

        return

    if status != "Completed": return

    # REFORMAT OUTPUT (wfcopy) AND GENERATE MANIFEST
    orig_dir = os.path.join(CROMWELL_EXECUTIONS_DIR, cromwell_id) # ?
    nice_dir = os.path.join()
    reformat_output(orig_dir, nice_dir):

    # TRANSFER RESULTS
    owner = user.User(user_id)
    transfer_client = owner.transfer_client()
    tdata = globus_sdk.TransferData(
        transfer_client,
        GLOBUS_ENDPOINT,
        dest_endpoint,
        label=run_id,
        sync_level="checksum",
        verify_checksum=True,
        preserve_timestamp=True,
        notify_on_succeeded=False, # TODO make an option?
        notify_on_failed=True,
        notify_on_inactive=True,
        skip_activation_check=False )
    tdata.add_item(nice_dir, dest_dir, recursive=True)
    transfer_result = transfer_client.submit_transfer(tdata)
    transfer_task_id = transfer_result["task_id"]

    # UPDATE DB RECORD
    # save task id


def reformat_output(cromwelldir, outdir):
    """
    Reformat cromwell run output and transfer to destination via Globus.
    """
    logdir = os.path.join(outdir, "log")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    cromwellFilesToSkip = ['stdout.background', 'stderr.background', 'script.background', 'script.submit']
    taskname = None
    rcfile = os.path.join(logdir, "workflow.rc")
    
    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, 'w') as fh:
        fh.write("#ExitCode\tTask\n")

    for rootdir, subdirs, files in os.walk(cromwelldir):
        if os.path.basename(rootdir).startswith('call-') and rootdir == os.path.join(cromwelldir, os.path.basename(rootdir)):
            taskname = re.sub(r'^call-', '', os.path.basename(rootdir))
            
        if rootdir.endswith('execution'):
            parentdir = str(pathlib.Path(rootdir).parent)
            shardname = None

            if re.search(r'shard-', parentdir):
                shardname = os.path.basename(parentdir)

            taskdir = os.path.join(outdir, taskname)
            if not os.path.exists(taskdir):
                os.makedirs(taskdir)

            for dname in subdirs:
                rsync(os.path.join(rootdir, dname), taskdir)

            for fname in files:
                fullname = os.path.join(rootdir, fname)
                outname = "%s-%s"%(taskname, shardname) if shardname else taskname

                if fname == 'stdout':
                    rsync(fullname, os.path.join(logdir, "%s.stdout"%outname))
                elif fname == 'stderr':
                    rsync(fullname, os.path.join(logdir, "%s.stderr"%outname))
                elif fname == 'script':
                    rsync(fullname, os.path.join(logdir, "%s.script"%outname))
                elif fname == 'rc':
                    with open(fullname, 'r') as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, 'a') as fh:
                        fh.write("%s\t%s\n"%(exitcode, outname))
                elif fname not in cromwellFilesToSkip:
                    rsync(fullname, taskdir))
                    
def runCommand(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.strip(), stderr.strip(), process.returncode

def rsync(src, dest):
    cmd = "rsync -a %s %s"%(src, dest)
    if DEBUG: print(cmd)
    stdout, stderr, exitcode = runCommand(cmd)
    if exitcode:
        sys.stderr.write(stderr)
        sys.exit(exitcode)





## SCHEDULE

# Examples:
#schedule.every(10).seconds.do(hello)
#schedule.every(1).minutes.do(process_completed_runs)
#schedule.every(1).hour.do(something_else)
#schedule.every(1).day.at("23:00").do(some_daily_task)
#schedule.every().sunday.at("4:00").do(some_weekly_task)


schedule.every(5).minutes.do(check_uploads)
schedule.every(5).minutes.do(check_runs)


while True:
    schedule.run_pending()
    time.sleep(10)
