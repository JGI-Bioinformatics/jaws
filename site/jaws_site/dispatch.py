#!/usr/bin/env python

"""
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import sys
import os
import json
import requests
from http.client import responses
from jaws_site import models

DEBUG = True if "JAWS_DEBUG" in os.environ else False
CROMWELL_URL = os.environ["CROMWELL_URL"]
WORKFLOWS_URL = "%s/api/workflows/v1" % (CROMWELL_URL,)
ENGINE_URL = "%s/engine/v1" % (CROMWELL_URL,)

class Dispatcher:
    """
    Analysis class for JAWS Site features function dispatch.  Each method must return JSON-RPC2 compliant response.
    """

    def dispatch(self, method, params):
        """
        Execute appropriate function, and return dict, which is either:
            response = { "result" : <string|integer|dict|list> }
        or
            response = { "error" : { "code" : <integer>, "message" : <string> }}
        where "message" is optional 
        """
        # TODO THIS SHOULD BE A DISPATCH TABLE
        if method == "server_status": return server_status(params)
        elif method == "task_status": return task_status(params)
        elif method == "run_metadata": return run_metadata(params)
        elif method == "get_labels": return get_labels(params)
        elif method == "cancel_run": return cancel_run(params)
        else:
            if DEBUG: print("Unknown method: %s" % (method,))
            return failure(400, "Unknown method")


def success(result): return { "jsonrpc" : "2.0", "result" : result }

def failure(code, message=None):
    if message is None:
        message = responses["status_code"] if "status_code" in responses else "Unknown error"
    return { "jsonrpc" : "2.0", "error" : { "code" : code, "message" : message }}


## METHODS


def server_status(params):
    """
    Return the current health status of any monitored subsystems
    """
    url = "%s/status" % (ENGINE_URL,)
    if DEBUG: print("Checking Cromwell status at %s" % (url,))
    try:
        r = requests.get(url)
    except:
        return failure(503, "Cromwell server timeout")
    if DEBUG: print("\tCromwell returned %s" % (r.status_code,))
    if r.status_code != requests.codes.ok: return failure(r.status_code)
    result = { "Cromwell" : "UP" }
    return success(result)


def run_metadata(params):
    """
    Retrieve the metadata of a run.    Neither authentication nor ownership is required.
    """
    if "cromwell_id" not in params: return failure(400)
    url = "%s/%s/metadata" % (WORKFLOWS_URL, params["cromwell_id"])
    try:
        r = requests.get(url)
    except:
        return failure(500, "Cromwell server timeout")
    if r.status_code != requests.codes.ok: return failure(r.status_code)
    return success(r.json())


def task_status(params):
    """
    Returns a list of task:status tuples, ordered by start time.
    """
    if "cromwell_id" not in params: return failure(400)
    url = "%s/%s/metadata" % (WORKFLOWS_URL, params["cromwell_id"])
    try:
        r = requests.get(url)
    except:
        return failure(500, "Cromwell server timeout")
    if r.status_code != requests.codes.ok: return failure(r.status_code)
    metadata = r.json()

    tasks = {}
    for task_name in metadata['calls']:
        start = ''
        try:
            start = metadata['calls'][task_name][0]['start']
        except:
            start = '?'
        end = ''
        try:
            end = metadata['calls'][task_name][0]['start']
        except:
            end = '?'
        tasks[start] = ( task_name, metadata['calls'][task_name][0]['executionStatus'], start, end )

    result = []
    for start in reversed(sorted(tasks.keys())): result.append(tasks[start])
    return success(result)


def cancel_run(params):
    """
    Cancel a run.
    """
    if "cromwell_id" not in params.keys(): return failure(400)
    url = "%s/%s/abort" & ( WORKFLOWS_URL, params["cromwell_id"] )
    try:
        r = requests.post(url)
    except:
        return failure(500, "Cromwell server timeout")
    if r.status_code != requests.codes.ok: return failure(r.status_code)
    result = r.json()
    result["status"] = "Cancelling"
    return success(result)


def get_labels(params):
    """
    Retrieve labels for a run.
    """
    if "cromwell_id" not in params: return failure(400)
    url = "%s/%s/labels" % (WORKFLOWS_URL, params["cromwell_id"])
    try:
        r = requests.get(url)
    except:
        return failure(500, "Cromwell server timeout")
    if r.status_code != requests.codes.ok: return failure(r.status_code)
    result = r.json()
    labels = result['labels']
    del labels['cromwell-workflow-id']
    return success(labels)
