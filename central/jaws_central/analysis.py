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
from time import sleep
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

import rpc_client
import models
from config import db

# DB MODELS
user_schema = models.UserSchema()
users_schema = models.UserSchema(many=True)
run_schema = models.RunSchema()
runs_schema = models.RunSchema(many=True)


# JAWS-SITE CONFIG
if "JAWS_SITES_CONFIG" not in os.environ: sys.exit('Env var $JAWS_SITES_CONFIG not defined')
site_config = configparser.ConfigParser()
site_config.read_file(open(os.environ["JAWS_SITES_CONFIG"]))

# INIT SITE RPC OBJECTS
site_rpc = {}    # name : rpc_client object
for site_name in site_config.sections():
    print("Initializing RPC client for %s" % (site_name,))
    site_rpc[site_name] = rpc_client.RPC_Client(
        site_config[site_name]["amqp_host"],
        site_config[site_name]["amqp_user"],
        site_config[site_name]["amqp_password"],
        site_config[site_name]["amqp_queue"] )


def _rpc_call(rpc_client, method, params={}):
    """
    Generic RPC using RPC-JSON2
    """
    # MAKE JSON-RPC2 PAYLOAD STRING
    query = {
        "jsonrpc" : "2.0",
        "method" : method,
        "params" : params
    }
    payload = json.dumps(query)

    # Send the request and store the requests Unique ID.
    corr_id = rpc_client.send_request(payload)

    # Wait until we have received a response.
    response = {}
    wait_interval = 0.1
    max_wait = 10
    waited = 0
    while rpc_client.queue[corr_id] is None:
        waited = waited + wait_interval    
        if waited > max_wait:
            response["error"] = {
                "code" : 500,
                "message": "Server timeout"
            }
            break
        sleep(wait_interval)

    # Return the response to the user (may be error).
    response_string = rpc_client.queue[corr_id]
    try:
        response = json.loads(response_string)
    except:
        response["error"] = {
            "code" : 500,
            "message" : "Invalid response"
        }
    return response


def _query_all_sites(command, params):
    """
    Submit the same command to every site and aggregate the responses.  Returns error if any site has an error.
    """
    any_errors = False
    results = {}
    for site, rpc_client in site_rpc.items():
        a_response = _rpc_call(rpc_client, command, params)
        if "error" in a_response: any_errors = True
        results[site] = a_response
        if "jsonrpc" in a_response: del(a_response["jsonrpc"])
    response = {
        "jsonrpc" : "2.0",
        "result" : results
    }
    if any_errors is True:
        response["error"] = {
            "code" : 500,
            "message" : "At least one subsystem returned an error"
        }
    return response


def pick_site(user, preferred_site=None):
    """
    If a preferred_site is specified, compare run's parameters to the Site's resources; if the run exceeds the Site's capabilities, return an alternate Site insead.  If no preferred site is suggested, then pick one based upon both capabilities and availability.
    """
    max_ram_mb = request.form.get("max_ram_mb")
    transfer_gb = request.form.get("transfer_gb")

    if preferred_site is None:
        # SELECT SITE
        # TODO NYI
        preferred_site = "lbnl"

    if preferred_site is not None:
        if max_ram_mb > site_config[preferred_site]["max_ram_mb"]: abort(400, "The required RAM exceeds the requested Site's available RAM")
        if transfer_gb > site_config[preferred_site]["max_transfer_gb"]: abort(400, "The input size exceeds the Site's data transfer limits")
        result = { "site" : preferred_site }
        return result

    return result


def site_info(user, site):
    """
    Returns Globus info for a JAWS-Site, so Client can send it a run's inputs.
    """
    if not site_config.has_section(site): abort(404, "Site not found")
    # TODO CHECK IF SITE IS ONLINE; OTHERWISE abort(503, "The requested resource is currently unavailable")

    # RETURN ONLY GLOBUS INFO, NOT AMQP SETTINGS
    result = {}
    result["endpoint"] = site_config[site]["globus_endpoint"]
    result["basepath"] = site_config[site]["globus_basepath"]
    result["staging_subdir"] = site_config[site]["staging_subdir"]
    return result


def _find_run(submission_id):
    """
    Returns the computing site name, run id, and owner user matching the submission id; aborts if not found.
    """
    run = models.Run.query.get(submission_id)
    if run is None: abort(404, "Run not found")
    return run


def _run_rpc_call(user, submission_id, method, params={}, restricted=True):
    """
    Generic RPC for analysis runs -- maps submission_id to a site and makes the RPC call using the corresponding RPC client object.  Uses JSON-RPC2.  Some methods are restricted to the owner while others may be executed by any valid user.
    """
    about = _find_run(submission_id)
    if "site" not in about:
        abort(404, "Run not found")
    if restricted is True:
        if about["owner"] == user: # TODO add admin pass
            pass
        else:
            abort(401, "Access denied")
    params["run_id"] = about["run_id"]
    site = about["site"]
    rpc_client = site_rpc[site]
    response = _rpc_call(rpc_client, method, params)
    if "error" in response:
        error = response["error"]
        abort(error["code"], error["message"])
    return response["result"], 200


def user_queue(user):
    """
    Return the current user's unfinished runs.
    """
    u = models.User.query.get(user)
    username = u.name
    return _query_all_sites("user_queue", {"username":username})


def search_runs(user):
    """
    Search for runs matching some criteria
    """
    params = {}
    query_user = request.form.get("query_user")
    if query_user:
        query_user = query_user.lower()
        if query_user != "all":
            params["uid"] = query_user
    else:
        params["uid"] = user
    status = request.form.get("status")
    if status:
        params["status"] = status
    delta_days = request.form.get("delta_days")
    if delta_days:
        d = datetime.today() - timedelta(int(delta_days))
        start_date = d.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        params["start"] = start_date
    start_date = request.form.get("start_date")
    if start_date:
        params["start"] = start_date
    end_date = request.form.get("end_date")
    if end_date:
        params["end"] = end_date
    result = {}
    site = request.form.get("site")
    if site:
        if site in site_rpc:
            rpc_client = site_rpc(site)
            result[site] = _rpc_call(rpc_client, "search_runs", params)
        else:
            abort(404, "Not a valid site")
    else:
        for site, rpc_client in site_rpc.items():
            result[site] = _rpc_call(rpc_client, "search_runs", params)
    return result


def submit_run(user, site, globus_transfer_task_id, submission_uuid):
    """
    Record the run submission in the RDB and relay to Site.
    Central stores minimal run metadata in it's database.
    Site will wait for the transfer to be complete before submitting to Cromwell.
    """
    if not site_config.has_section(site): abort(404, "Site not found")

    # SAVE MINIMAL RUN METADATA IN CENTRAL DB
    # SUBMISSION UUID AND TRANSFER_TASK_ID ARE NOT STORED;
    # THEY ARE STORED BY THE SITE INSTEAD
    submission_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    stmt = "INSERT INTO submissions (owner, site, submission_date, state) VALUES (%s,%s,%s, 'INITIATING')"

    # TODO CHANGE TO SQLALCHEMY
    cnx = current_app.get_db()
    cursor = cnx.cursor()
    try:
        cursor.execute(stmt, (user, site, submission_date))
    except Error as e:
        abort(500, "Error inserting into db")
    run_id = cursor.lastrowid

    # GET USER'S GLOBUS AUTH TOKEN
    auth_token = "TODO" # TODO

    # SUBMIT TO SITE
    params = {
        "auth_token" : auth_token,
        "run_id" : run_id,
        "submission_uuid" : submission_uuid,
        "transfer_task_id" : globus_transfer_task_id
    }
    rpc_client = site_rpc[site]
    response = _rpc_call(rpc_client, method, params)
    if "error" in response:
        error = response["error"]
        # CHANGE RUN STATE TO "SUBMISSION_FAILED"
        # TODO
        abort(error["code"], error["message"])

    result = response["result"]
    result["run_id"] = run_id
    return result, 201


def run_status(user, submission_id):
    """
    Retrieve the current status of a run.
    """
    return _run_rpc_call(user, submission_id, "run_status", restricted=False)


def task_status(user, submission):
    """
    Retrieve run status with task-level detail.
    """
    return _run_rpc_call(user, submission_id, "task_status", restricted=False)


def run_metadata(user, submission_id):
    """
    Retrieve the metadata of a run.
    """
    return _run_rpc_call(user, submission_id, "run_metadata", restricted=False)


def get_labels(submission_id):
    """
    Retrieve labels for a run.
    """
    return _run_rpc_call(user, submission_id, "get_labels", restricted=False)


def cancel_run(user, submission_id):
    """
    Cancel a run.
    """
    return _run_rpc_call(user, submission_id, "cancel", params={}, restricted=True)


def invalidate_task(user, submission_id, task_name):
    """
    Mark a task's output as invalid to prevent it from being cached.
    """
    params = { "task" : task_name }
    return _run_rpc_call(user, submission_id, "invalidate", params, restricted=True)
