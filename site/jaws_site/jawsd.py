#!/usr/bin/env python3

"""
JAWS Daemon process periodically checks on runs and performs actions to usher them to the next state.
"""

import sys
import schedule
import time
import os
import requests
import json
import globus_sdk
#from globus_sdk import AuthClient, AccessTokenAuthorizer

from jaws_site import wfcopy, models


DEBUG = True if "JAWS_DEBUG" in os.environ else False

# SITE
STAGING_DIR = os.environ["JAWS_STAGING_DIR"]
RESULTS_DIR = os.environ["JAWS_RESULTS_DIR"]

# CROMWELL
CROMWELL_URL = os.environ["CROMWELL_URL"]
WORKFLOWS_URL = "%s/api/workflows/v1" % (CROMWELL_URL,)

# GLOBUS
GLOBUS_CLIENT_ID = os.environ["GLOBUS_CLIENT_ID"]
AUTH_SERVICE_NAME = 'auth.globus.org'
TRANSFER_SERVICE_NAME = 'transfer.api.globus.org'
GROUPS_SERVICE_NAME = "04896e9e-b98e-437e-becd-8084b9e234a0" # NOTE NOT groups.api.globus.org !?!
globus_client = globus_sdk.NativeAppAuthClient(GLOBUS_CLIENT_ID)


def globus_transfer_client(user):
    """
    Get Globus transfer client object for user.
    """
    access_token = user.transfer_access_token
    authorizer = globus_sdk.AccessTokenAuthorizer(access_token)
    transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
    return transfer_client


def check_uploads():
    """
    Check Globus transfer status for uploading runs.
    """
    q = db.session.query(Run).filter_by(status="uploading")
    for run in q:
        check_transfer_status(run)


def check_transfer_status(run):
    """
    Check on the status of one uploading run.
    """
    if DEBUG: print("Uploading: %s" % (run.id,))
    user = db.session.query(User).get(run.user_id)
    transfer_client = globus_transfer_client(user)
    task = transfer_client.get_task(transfer_task_id)
    globus_status = task["status"]
    if globus_status == "ACTIVE":
        pass
    elif globus_status == "FAILED":
        run.status = "upload failed"
        db.session.commit()
    elif globus_status == "INACTIVE":
        run.status = "upload inactive"
        db.session.commit()
    elif globus_status == "SUCCEEDED":
        submit_run(run)


def submit_run(run):
    """
    Submit a run to Cromwell.
    """
    # DEFINE FILES
    wdl_file = os.path.join(STAGING_DIR, run.submission_uuid+".wdl")
    json_file = os.path.join(STAGING_DIR, run.submission_uuid+".json")
    zip_file = os.path.join(STAGING_DIR, run.submission_uuid+".zip") # may not exist

    # POST TO CROMWELL SERVER
    files = {}
    files['workflowInputs'] = ( 'workflowInputs', open(json_file, 'r'), 'application/json' )
    files['workflowSource'] = ( 'workflowSource', open(wdl_file, 'r'), 'application/json' )
    if os.path.exists(zip_file):
        files['workflowDependencies'] = ( 'workflowDependencies', open(zip_file, 'rb'), 'application/zip' )
    try:
        r = requests.post(WORKFLOWS_URL, files=files)
    except:
        run.status = "cromwell timeout"
        db.session.commit()
        return
    if r.status_code == 201:
        run.cromwell_id = result.json()[0]["id"]
        run.status = "submitted"
    else:
        run.status = "submission failed"
    db.session.commit()


def check_cromwell():
    """
    Check current runs that have been submitted to Cromwell.
    """
    q = db.session.query(Run).filter(Run.status.in_(["submitted", "queued", "running"]))
    for run in q:
        check_run_status(run)


def check_run_status(run):
    """
    Check on the status of one Run.
    """
    if DEBUG: print("Queued/Running (Cromwell): %s" % (run.id,))
    try:
        r = requests.get("%s/%s/status" % (CROMWELL_URL, run.cromwell_id))
    except:
        # Cromwell server timeout
        return
    if r.status_code != requests.codes.ok: return

    cromwell_status = r.json()["status"]
    if cromwell_status == "Failed":
        run.status = "failed"
        db.session.commit()
    elif cromwell_status == "Completed":
        finish_run(run)


def finish_run(run):
    """
    Reformat and send run output.
    """
    # GET OUTPUT PATH
    url = "%s/%s/metadata" % (WORKFLOWS_URL, run.cromwell_id)
    try:
        r = requests.get(url)
    except:
        # Cromwell server timeout
        return
    if r.status_code != requests.codes.ok: return
    orig_dir = r.json()["workflowRoot"]

    # UPDATE STATE, RUN WFCOPY
    run.status = "post-processing"
    db.session.commit()
    nice_dir = os.path.join(RESULTS_DIR, run.id)
    wfcopy.wfcopy(orig_dir, nice_dir, flattenShardDir=False, verbose=False)

    # UPDATE STATE, SUBMIT XFER TO GLOBUS
    run.status = "downloading output"
    db.session.commit()
    user = db.session.query(User).get(run.user_id)
    transfer_client = globus_transfer_client(user)
    tdata = globus_sdk.TransferData(
        transfer_client,
        GLOBUS_ENDPOINT,
        run.dest_endpoint,
        label=run.id,
        sync_level="checksum",
        verify_checksum=True,
        preserve_timestamp=True,
        notify_on_succeeded=True, # TODO make an option?
        notify_on_failed=True,
        notify_on_inactive=True,
        skip_activation_check=False )
    tdata.add_item(nice_dir, run.dest_path, recursive=True)
    transfer_result = transfer_client.submit_transfer(tdata)
    run.download_task_id = transfer_result["task_id"]
    db.session.commit()


if __name__ == "__main__":
    """
    The daemon is usually started by server.py, but can be executed directly for debugging purposes.
    """
    # NOTE: the real schedule is specified in server.py; this is for debugging only:
    schedule.every(5).seconds.do(check_uploads)
    schedule.every(5).seconds.do(check_cromwell)
    while True:
        schedule.run_pending()
        time.sleep(1)
