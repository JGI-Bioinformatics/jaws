#!/usr/bin/env python

"""
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
from jaws_site import config


class Dispatcher:
    """
    Analysis class for JAWS Site features function dispatch.  Each method must return JSON-RPC2 compliant response.
    """

    cromwell_engine_status_url = None
    workflows_url = None

    def __init__(self):
        """
        Init obj
        """
        self.logger = logging.getLogger(__package__)
        conf = config.JawsConfig()
        self.cromwell_engine_status_url = conf.get_cromwell('engine_status_url')
        self.workflows_url = conf.get_cromwell("workflows_url")

    def dispatch(self, method, params):
        """
        Execute appropriate function, and return dict, which is either:
            response = { "result" : <string|integer|dict|list> }
        or
            response = { "error" : { "code" : <integer>, "message" : <string> }}
        where "message" is optional
        """
        # TODO THIS SHOULD BE A DISPATCH TABLE
        if method == "server_status":
            return self.server_status(params)
        elif method == "task_status":
            return self.task_status(params)
        elif method == "run_metadata":
            return self.run_metadata(params)
        elif method == "get_labels":
            return self.get_labels(params)
        elif method == "cancel_run":
            return self.cancel_run(params)
        else:
            self.logger.warning(f"Invalid RPC method: {method}")
            return self.failure(400, "Unknown method")

    def success(self, result):
        """
        Return a JSON-RPC2 successful result message.
        """
        return {"jsonrpc": "2.0", "result": result}

    def failure(self, code, message=None):
        """
        Return a JSON-RPC2 error message.
        """
        if message is None:
            message = (
                responses["status_code"]
                if "status_code" in responses
                else "Unknown error"
            )
        return {"jsonrpc": "2.0", "error": {"code": code, "message": message}}

    def server_status(self, params):
        """
        Return the current health status of any monitored subsystems
        """
        try:
            r = requests.get(self.cromwell_engine_status_url)
        except Exception:
            self.logger.warning("Cromwell server timeout")
            return self.failure(503, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        result = {"Cromwell": "UP"}
        return self.success(result)

    def run_metadata(self, params):
        """
        Retrieve the metadata of a run.    Neither authentication nor ownership is required.
        """
        if "cromwell_id" not in params:
            self.logger.error(f"run_metadata rpc did not include cromwell_id")
            return self.failure(400)
        url = "%s/%s/metadata" % (self.workflows_url, params["cromwell_id"])
        try:
            r = requests.get(url)
        except Exception:
            self.logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        return self.success(r.json())

    def task_status(self, params):
        """
        Returns a list of task:status tuples, ordered by start time.
        """
        if "cromwell_id" not in params:
            self.logger.error(f"task_status rpc did not include cromwell_id")
            return self.failure(400)
        url = "%s/%s/metadata" % (self.workflows_url, params["cromwell_id"])
        try:
            r = requests.get(url)
        except Exception:
            self.logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        metadata = r.json()

        tasks = {}
        for task_name in metadata["calls"]:
            start = ""
            try:
                start = metadata["calls"][task_name][0]["start"]
            except Exception:
                start = "?"
            end = ""
            try:
                end = metadata["calls"][task_name][0]["start"]
            except Exception:
                end = "?"
            tasks[start] = (
                task_name,
                metadata["calls"][task_name][0]["executionStatus"],
                start,
                end,
            )

        result = []
        for start in reversed(sorted(tasks.keys())):
            result.append(tasks[start])
        return self.success(result)

    def task_ids(self, params):
        """
        Returns a list of all JTM task ids associated with a run.  May not be complete if run in progress.
        If a task was rerun automatically by Cromwell, only the last attempt is included (i.e. no redundant tasks).
        """
        if "cromwell_id" not in params:
            self.logger.error(f"task_ids rpc did not include cromwell_id")
            return self.failure(400)
        url = "%s/%s/metadata" % (self.workflows_url, params["cromwell_id"])
        try:
            r = requests.get(url)
        except Exception:
            self.logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        metadata = r.json()
        tasks = []
        for task_name in metadata["calls"]:
            try:
                jobId = metadata["calls"][task_name][-1]["jobId"]
                tasks.append(jobId)
            except Exception:
                pass
        return self.success(tasks)

    def cancel_run(self, params):
        """
        Cancel a run.
        """
        if "cromwell_id" not in params.keys():
            self.logger.error(f"cancel_run rpc did not include cromwell_id")
            return self.failure(400)
        url = "%s/%s/abort" & (self.workflows_url, params["cromwell_id"])
        try:
            r = requests.post(url)
        except Exception:
            self.logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        result = r.json()
        result["status"] = "Cancelling"
        return self.success(result)
