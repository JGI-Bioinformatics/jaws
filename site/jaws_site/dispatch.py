"""
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
from jaws_site import config


class Dispatcher:
    """
    Analysis class for JAWS Site features function dispatch.
    Each method must return JSON-RPC2 compliant response.
    """

    cromwell_engine_status_url = None
    workflows_url = None

    def __init__(self):
        """Constructor"""
        self.logger = logging.getLogger(__package__)
        conf = config.JawsConfig()
        self.cromwell_engine_status_url = conf.get_cromwell('engine_status_url')
        self.workflows_url = conf.get_cromwell("workflows_url")

    def dispatch(self, method, params):
        """Given a method keyword, call the apporpriate function with the provided parameters.

        :param method: The method requested.  Abort if unrecognized.
        :type method: str
        :param params: Associated parameters, if any.  Varies by method.
        :type params: dict
        :return: Response in JSON-RPC2 format.
        :rtype: dict
        """
        # TODO THIS SHOULD BE A DISPATCH TABLE
        if method == "server_status":
            return self.server_status()
        elif method == "task_status":
            return self.task_status(params.get("cromwell_id", default=None))
        elif method == "run_metadata":
            return self.run_metadata(params.get("cromwell_id", default=None))
        elif method == "get_labels":
            return self.get_labels(params.get("cromwell_id", default=None))
        elif method == "cancel_run":
            return self.cancel_run(params.get("cromwell_id", default=None))
        else:
            self.logger.warning(f"Invalid RPC method: {method}")
            return self.failure(400, "Unknown method")

    def success(self, result):
        """Return a JSON-RPC2 successful result message.

        :param result: The result returned by a successful RPC call.
        :type result: str|int|dict|list
        :return: JSON-RPC2 formatted response
        :rtype: str
        """
        return {"jsonrpc": "2.0", "result": result}

    def failure(self, code, message=None):
        """Return a JSON-RPC2 error message.

        :param code: The error code returned by the RPC call
        :type code: int
        :param message: The error message returned by the RPC call
        :type message: str, optional
        :return: JSON-RPC2 formatted response
        :rtype: dict
        """
        if message is None:
            message = (
                responses["status_code"]
                if "status_code" in responses
                else "Unknown error"
            )
        return {"jsonrpc": "2.0", "error": {"code": code, "message": message}}

    def server_status(self):
        """Return the current status of the Cromwell server.

        :return: Either a success- or failure-formatted JSON-RPC2 response, if Cromwell up or not.
        :rtype: dict
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

    def run_metadata(self, cromwell_id):
        """Retrieve the metadata of a run.

        :param cromwell_id: Cromwell run ID
        :type params: str
        :return: The Cromwell metadata for the specified run.
        :rtype: dict
        """
        if not cromwell_id:
            self.logger.error(f"run_metadata rpc did not include cromwell_id")
            return self.failure(400)
        url = f'{self.workflows_url}/{cromwell_id}/metadata'
        try:
            r = requests.get(url)
        except Exception:
            self.logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        return self.success(r.json())

    def task_status(self, cromwell_id):
        """Returns the status for each task processed thus far in a run.

        :param cromwell_id: The Cromwell run ID
        :type cromwell_id: str
        :return: A list of task:status tuples, ordered by start time.
        :rtype: list
        """
        if not cromwell_id:
            self.logger.error(f"task_status rpc did not include cromwell_id")
            return self.failure(400)
        url = f'{self.workflows_url}/{cromwell_id}/metadata'
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

    def task_ids(self, cromwell_id):
        """Returns a list of all JTM task ids associated with a run thus far.

        :param cromwell_id: The Cromwell run ID
        :type cromwell_id: str
        :return: A list of JTM task IDs.  If a task was rerun, only the last is included.
        :rtype: list
        """
        if not cromwell_id:
            self.logger.error(f"task_ids rpc did not include cromwell_id")
            return self.failure(400)
        url = f'{self.workflows_url}/{cromwell_id}/metadata'
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

    def cancel_run(self, cromwell_id):
        """Cancel a run.

        :param cromwell_id: The Cromwell run ID
        :type cromwell_id: str
        :return: Either a JSON-RPC2-compliant success or failure message, if the run could be cancelled or not.
        :rtype: dict
        """
        if not cromwell_id:
            self.logger.error(f"cancel_run rpc did not include cromwell_id")
            return self.failure(400)
        url = f'{self.workflows_url}/{cromwell_id}/abort'
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
