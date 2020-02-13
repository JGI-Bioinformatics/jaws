#!/usr/bin/env python

"""
Run management.  This module provides REST endpoints and relays/aggregates communications with the JAWS Sites via AMQP-RPC.  RPC communications use JSON-RPC2.
"""

from datetime import datetime, timedelta
import sys
import os
from flask import make_response, abort, Flask, request, redirect, url_for
from flask import current_app
#from werkzeug.utils import secure_filename
import time
#import glob
#import subprocess
#from subprocess import Popen, PIPE
#import markdown
import json
import requests
from requests import Session
#import tempfile
#import shutil
import mysql.connector
from mysql.connector import Error
#import socket
#import re
import threading
import amqpstorm
from amqpstorm import Message
import uuid
from datetime import datetime
import configparser
import csv

import rpc_client
import models
from config import db
#import user as jaws_user

from models import Site, Run, User, Workflow

DEBUG = True if "JAWS_DEBUG" in os.environ else False


# JAWS-SITE CONFIG
if "JAWS_SITES_CONFIG" not in os.environ: sys.exit('Env var $JAWS_SITES_CONFIG not defined')
site_config = configparser.ConfigParser()
site_config.read_file(open(os.environ["JAWS_SITES_CONFIG"]))


# INIT SITE RPC OBJECTS
rpc_clients = {} # site name : rpc client object
for site_name in site_config.sections():
    print("Initializing RPC client for %s" % (site_name,))
    rpc_clients[site_name] = rpc_client.RPC_Client(
        site_config[site_name]["amqp_host"],
        site_config[site_name]["amqp_user"],
        site_config[site_name]["amqp_password"],
        site_config[site_name]["amqp_queue"],
        vhost=site_config[site_name]["amqp_vhost"] )



## DB FUNCTIONS -- faster/less overhead than RPC, particularly when multiple sites are involved in the query (e.g. user_queue)

def user_queue(user):
    """
    Return the current user's unfinished runs.
    """
    result = db.session.query(Run).filter_by(user_id=user).filter(Run.status.in_(["uploading","queued","running","post-processing","downloading"])).all()
    return result


def user_recent(user):
    """
    Return the current user's recent runs, regardless of status.
    """
    delta_days = request.form.get("delta_days", 30)
    start_date = datetime.today() - timedelta(int(delta_days))
    result = db.session.query(Run).filter_by(user_id=user).filter(Run.submitted >= start_date).all()
    return result


def get_site(user):
    """
    If a preferred_site is specified, compare run's parameters to the Site's resources; if the run exceeds the Site's capabilities, return an alternate Site insead.  If no preferred site is suggested, then pick one based upon both capabilities and availability.
    """
    # TODO PICK SITE
    # TODO COMPARE REQUIRED RESOURCES TO SITE AVAILABLE RESOURCES
    # TODO CHECK IF SITE IS ONLINE; OTHERWISE abort(503, "The requested resource is currently unavailable")

    site_id = request.form.get("site", None)
    if site_id: site_id = site_id.upper()
    else: site_id = "LBNL"

    max_ram_gb = request.form.get("max_ram_gb", None)
    transfer_gb = request.form.get("transfer_gb", None)

    if DEBUG: print("Return site info for: %s" % (site_id,))
    site = db.session.query(Site).get(site_id)
    if not site: abort(404, "Site not found")

    # RETURN GLOBUS INFO
    result = { "endpoint" : site.endpoint, "site" : site.id, "staging" : site.staging }
    #result["dir"] = os.path.join(site.basepath, site.staging)
    return result, 200


def submit_run(user):
    """
    Record the run submission in the database.  The status of the new run will be "uploading".
    """
    if DEBUG: print("Init run submission from %s" % (user,))

    # GET PARAMS
    site_id = request.form.get("site", None).upper()
    submission_uuid = request.form.get("submission_uuid", None)
    globus_transfer_task_id = request.form.get("transfer_task_id", None)
    globus_endpoint = request.form.get("globus_endpoint", None)
    outdir = request.form.get("outdir", None)

    # VALIDATE
    site = db.session.query(Site).get(site_id)
    if not site: abort(404, "Site not found")
    # TODO VERIFY GLOBUS TRANSFER HAS NOT FAILED

    # SAVE IN DB
    if DEBUG: print("Inserting into database")
    current_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    run = Run(
        user_id=user,
        site_id=site_id,
        status="uploading",
        submission_uuid=submission_uuid,
        upload_task_id=globus_transfer_task_id,
        dest_endpoint = globus_endpoint,
        dest_path = outdir )
    db.session.add(run)
    db.session.commit()

    # RETURN PRIMARY KEY
    result = { "run_id" : run.id }
    return result, 201


def run_status(user, run_id):
    """
    Retrieve the current status of a run.
    """
    q = db.session.query(Run).get(run_id)
    if not q: abort(404, "Run not found")
    result = { "status" : q.status }
    return result, 200



## RPC FUNCTIONS

def _rpc(user, run_id, method, params = {}):
    """
    Performs the specified RPC function and returns result if OK, aborts if error.
    """
    run = db.session.query(Run).get(run_id)
    if not run: abort(404, "Run not found")
    if run.user_id != user: abort(401, "Access denied")
    rpc = rpc_clients[run.site_id] 
    #if run.cromwell_id is None: abort(400, "Run not yet scheduled")
    params["cromwell_id"] = run.cromwell_id
    result = rpc.request(method, params)
    if "error" in result: abort(result["error"]["code"], result["error"]["message"])
    return result, 200


def task_status(user, run_id):
    """
    Retrieve run status with task-level detail.
    """
    return _rpc(user, run_id, "task_status")


def run_metadata(user, run_id):
    """
    Retrieve the metadata of a run.
    """
    return _rpc(user, run_id, "run_metadata")


def cancel_run(user, run_id):
    """
    Cancel a run.
    """
    return _rpc(user, run_id, "cancel_run")


#def invalidate_task(user, run_id, task_name):
#    """
#    Mark a task's output as invalid to prevent it from being cached.
#    """
#    return _rpc(user, run_id, "invalidate_task", { "task" : task_name })


def get_labels(run_id):
    """
    Retrieve labels for a run.
    """
    return _rpc(user, run_id, "get_labels")


#def search_runs(user):
#    """
#    Search for runs matching some criteria
#    """
#    params = {}
#    query_user = request.form.get("query_user", None)
#    if query_user:
#        query_user = query_user.lower()
#        if query_user != "all":
#            params["uid"] = query_user
#    else:
#        params["uid"] = user
#    status = request.form.get("status", None)
#    if status:
#        params["status"] = status
#    delta_days = request.form.get("delta_days", 90)
#    if delta_days:
#        d = datetime.today() - timedelta(int(delta_days))
#        start_date = d.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
#        params["start"] = start_date
#    start_date = request.form.get("start_date", None)
#    if start_date:
#        params["start"] = start_date
#    end_date = request.form.get("end_date", None)
#    if end_date:
#        params["end"] = end_date
#    result = {}
#    site = request.form.get("site", None)
#    if site:
#        if site in site_rpc:
#            rpc_client = site_rpc(site)
#            result[site] = _rpc_call(rpc_client, "search_runs", params)
#        else:
#            abort(404, "Not a valid site")
#    else:
#        for site, rpc_client in site_rpc.items():
#            result[site] = _rpc_call(rpc_client, "search_runs", params)
#    return result
