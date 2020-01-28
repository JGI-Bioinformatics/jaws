#!/usr/bin/env python

"""
RPC functions for analysis/run management and reporting.  Most of these simple wrap the localhost Cromwell REST functions.
"""

import sys
import os
import time
from datetime import datetime, timedelta
import json
import requests
import tempfile
import mysql.connector
from http.client import responses

from jaws_site import config


class Analysis:
    """
    Analysis class for JAWS Site features function dispatch.
    """

    def __init__(self, config):
        """
        Constructor
        """
        print("Initialize analysis dispatcher") # ECCE
        self.config = config
        self.cromwell = "%s/api/workflows/v1" % (config["CROMWELL"]["url"],)
        self.engine = "%s/api/engine/v1" % (config["CROMWELL"]["url"],)
#        self.dispatch = {
#            "server_status" : self.server_status,
#            "user_queue" : self.user_queue
#        }
        
    def dispatch(self, body):
        """
        Unpack message body, execute appropriate function, and return result.
        """
        request = json.loads(body)
        method = request["method"]
        params = request["params"]

        if method == "server_status":
            return self.server_status(params)
        elif method == "user_queue":
            return self.user_queue(params)
        else:
            return self.failure(400, "Invalid request")

#    def __call__(self, method, params):
#        """
#        Dispatch method
#        """
#        return self.dispatch[method](params)

    def success(self, result):
        """
        Format responses according to JSON-RPC2 specification and return string.
        """
        response = {
            "jsonrpc" : "2.0",
            "result" : result
        }
        return json.dumps(response)

    def failure(self, status_code, message=None):
        """
        Format responses according to JSON-RPC2 specification and return string.
        """
        if not message:
            if status_code in responses:
                message = responses[status_code]
            else:
                message = "Unknown Error"
        response = {
            "jsonrpc" : "2.0",
            "error" : {
                "code" : status_code,
                "message" : message
            }
        }
        return json.dumps(response)


    def server_status(self, params):
        """
        Return the current health status of any monitored subsystems
        """
        r = requests.get("%s/status" % (self.engine,))
        if r.status_code == requests.codes.ok:
            return self.success(r.json())
        else:
            return self.failure(r.status_code)

    def user_queue(self, params):
        """
        Return a user's unfinished runs.
        """
        if "username" not in params: return self.failure(400)
        data = {
            "status" : [ 'Submitted', 'Running', 'Aborting'],
            "label" : [ "username:" + params["username"] ]
        }
        r = requests.get(self.cromwell+"/query", data)
        if r.status_code == requests.codes.ok:
            return self.success(r.json())
        else:
            return self.failure(r.status_code)

    def user_history(self, params):
        """
        Return a user's recent runs.
        """
        if not params.keys() >= {"delta_days","username"}: return self.failure(400)
        d = datetime.today() - timedelta(int(params["delta_days"]))
        start_date = d.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        data = {
            "label" : [ "username:" + params["username"] ],
            "start" : start_date
        }
        r = requests.get(self.cromwell+"/query", data)
        if r.status_code == requests.codes.ok:
            return self.success(r.json())
        else:
            return self.failure(r.status_code)

    def run_status(self, params):
        """
        Retrieve the current status of a run.
        """
        if "run_id" not in params: return self.failure(400)
        r = requests.get("%s/%s/status" % (self.cromwell, params["run_id"]))
        if r.status_code == requests.codes.ok:
            return self.success(r.json())
        else:
            return self.failure(r.status_code)

    def run_metadata(self, params):
        """
        Retrieve the metadata of a run.    Neither authentication nor ownership is required.
        """
        if "run_id" not in params: return self.failure(400)
        r = requests.get("%s/run_id/metadata" % (self.cromwell, params["run_id"]))
        if r.status_code == requests.codes.ok:
            return self.success(r.json())
        else:
            return self.failure(r.status_code)

    def task_status(self, params):
        """
        Returns a list of task:status tuples, ordered by start time.
        """
        if "run_id" not in params: return self.failure(400)
        r = requests.get("%s/run_id/metadata" % (self.cromwell, params["run_id"]))
        if not r.status_code == requests.codes.ok:
            return self.failure(r.status_code)
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
        for start in reversed(sorted(tasks.keys())):
            result.append(tasks[start])
        return self.success(result)

    def cancel_run(self, params):
        """
        Cancel a run.
        """
        if "run_id" not in params.keys(): return self.failure(400)
        r = requests.post(self.cromwell + "/run_id/abort", params["run_id"])
        if r.status_code == requests.codes.ok:
            result = r.json()
            result["status"] = "Cancelling"
            return self.success(result)
        else:
            return self.failure(r.status_code)

    def get_labels(self, params):
        """
        Retrieve labels for a run.
        """
        if "run_id" not in params: return self.failure(400)
        r = requests.get("%s/run_id/labels" % (self.cromwell, params["run_id"]))
        if r.status_code == requests.codes.ok:
            result = r.json()
            labels = result['labels']
            del labels['cromwell-workflow-id']
            return self.success(labels)
        else:
            return self.failure(r.status_code)

# TODO UPDATE LABELS
    
# TODO SEARCH

    # TODO ADD SUPPORT FOR LABELS YAML
    def submit_run(self, params):
        """
        Submit a runs.  Requires:
        - user uid
        - globus transfer task id
        - (local) paths to: main WDL, inputs JSON, and optionally subworkflows ZIP and labels YAML files
        """
        if not params.keys() >= {"inputs", "labels", "globus", "paths"}:
            return self.failure(400)

        # LABELS
        #files['labels'] = ('labels', open(labels_path, 'r'), 'application/yaml')
        tmp_labels_file = tmp_basename + ".labels.json"
        _prepare_labels_file(tmp_labels_file, username, labels)

        # DATA
        data = {}

        # SUBMIT TO CROMWELL
        files = {}
        files['workflowInputs'] = ( 'workflowInputs', open(json_file, 'r'), 'application/json' )
        files['workflowSource'] = ( 'workflowSource', open(wdl_file, 'r'), 'application/json' )
        files['labels'] = ( 'labels', open(tmp_labels_file, 'r'), 'application/json' )
        if zip_file is not None:
            files['workflowDependencies'] = ( 'workflowDependencies', open(zip_file, 'rb'), 'application/zip' )

        # SUBMIT RUN TO CROMWELL SERVER VIA REST
        r = requests.post(self.cromwell, data=data, files=files)
        if r.status_code == requests.codes.ok:
            return self.success(r.json())
        else:
            return self.failure(r.status_code)


def _prepare_labels_file(outfile, username, infile=None):
    """
    Given a username and optionally other labels (in a dict), create a JSON file.
    """
    assert outfile
    assert username
    labels = {}
    if infile is not None:
        with open(infile, 'r') as f:
            data = json.load(f)
            for key in data:
                labels[key] = data[key]
    labels["username"] = username
    with open(outfile, 'w') as f:
        json.dump(labels, f)

