"""
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
from typing import Dict
import os
import click
from jaws_site import config


logger = logging.getLogger(__package__)


class Dispatcher:
    """
    Analysis class for JAWS Site features function dispatch.
    Each method must return JSON-RPC2 compliant response.
    """

    cromwell_engine_status_url = None
    workflows_url = None

    def __init__(self, conf):
        """Constructor"""
        self.cromwell_engine_status_url = conf.get("CROMWELL", "engine_status_url")
        self.workflows_url = conf.get("CROMWELL", "workflows_url")

    def dispatch(self, method: str, params: Dict[str, str]):
        """Given a method keyword, call the apporpriate function with the provided parameters.

        :param method: The method requested.  Abort if unrecognized.
        :type method: str
        :param params: Associated parameters, if any.  Varies by method.
        :type params: dict
        :return: Response in JSON-RPC2 format.
        :rtype: dict
        """
        operations = {
            "server_status": self.server_status,
            "task_status": self.task_status,
            "run_metadata": self.run_metadata,
            "cancel_run": self.cancel_run,
            "run_logs": self.run_logs,
            "failure_logs": self.failure_logs,
        }
        proc = operations.get(method)
        if proc:
            return proc(params)
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

    def server_status(self, params):
        """Return the current status of the Cromwell server.

        :return: Either a success- or failure-formatted JSON-RPC2 response, if Cromwell up or not.
        :rtype: dict
        """
        try:
            r = requests.get(self.cromwell_engine_status_url)
        except requests.exceptions.RequestException:
            logger.warning("Cromwell server timeout")
            return self.failure(503, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        result = {"Cromwell": "UP"}
        return self.success(result)

    def run_metadata(self, params):
        """Retrieve the metadata of a run.

        :param cromwell_id: Cromwell run ID
        :type params: str
        :return: The Cromwell metadata for the specified run.
        :rtype: dict
        """
        if "cromwell_id" not in params:
            logger.error(f"run_metadata rpc did not include cromwell_id")
            return self.failure(400)
        url = f"{self.workflows_url}/{params['cromwell_id']}/metadata"
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        return self.success(r.json())

    def task_status(self, params):
        """Returns the status for each task processed thus far in a run.

        :param cromwell_id: The Cromwell run ID
        :type cromwell_id: str
        :return: A list of task:status tuples, ordered by start time.
        :rtype: list
        """
        if "cromwell_id" not in params:
            logger.error(f"task_status rpc did not include cromwell_id")
            return self.failure(400)
        url = f"{self.workflows_url}/{params['cromwell_id']}/metadata"
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        metadata = r.json()

        tasks = {}
        for task_name in metadata["calls"]:
            start = ""
            if "start" in metadata["calls"][task_name][0]:
                start = metadata["calls"][task_name][0]["start"]
            else:
                start = "?"
            end = ""
            if "end" in metadata["calls"][task_name][0]["start"]:
                end = metadata["calls"][task_name][0]["end"]
            else:
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
        """Returns a list of all JTM task ids associated with a run thus far.

        :param cromwell_id: The Cromwell run ID
        :type cromwell_id: str
        :return: A list of JTM task IDs.  If a task was rerun, only the last is included.
        :rtype: list
        """
        if "cromwell_id" not in params:
            logger.error(f"task_ids rpc did not include cromwell_id")
            return self.failure(400)
        url = f"{self.workflows_url}/{params['cromwell_id']}/metadata"
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        metadata = r.json()
        tasks = []
        for task_name in metadata["calls"]:
            if "jobId" in metadata["calls"][task_name][-1]:
                jobId = metadata["calls"][task_name][-1]["jobId"]
                tasks.append(jobId)
        return self.success(tasks)

    def cancel_run(self, params):
        """Cancel a run.

        :param cromwell_id: The Cromwell run ID
        :type cromwell_id: str
        :return: Either a JSON-RPC2-compliant success or failure message, if the run could be cancelled or not.
        :rtype: dict
        """
        if "cromwell_id" not in params:
            logger.error(f"cancel_run rpc did not include cromwell_id")
            return self.failure(400)
        url = f"{self.workflows_url}/{params['cromwell_id']}/abort"
        try:
            r = requests.post(url)
        except requests.exceptions.RequestException:
            logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        result = r.json()
        result["status"] = "Cancelling"
        return self.success(result)

    def run_logs(self, params):
        """
        Retrieve the Cromwell logs for a run.
        """
        if "cromwell_id" not in params:
            logger.error(f"run_logs rpc did not include cromwell_id")
            return self.failure(400)
        url = f'{self.workflows_url}/{params["cromwell_id"]}/logs'
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        return self.success(r.json())

    def failure_logs(self, params):
        """
        Concatenate stdout, stderr files of any failed tasks.
        """
        if "cromwell_id" not in params:
            self.logger.error(f"failure_logs rpc did not include cromwell_id")
            return self.failure(400)

        # GET RUN FOLDER
        url = f'{self.workflows_url}/{params["cromwell_id"]}/metadata'
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            self.logger.warning("Cromwell server timeout")
            return self.failure(500, "Cromwell server timeout")
        if r.status_code != requests.codes.ok:
            return self.failure(r.status_code)
        run_dir = r.json()["workflowRoot"]

        # FIND FAILED TASKS AND CONCATENATE TARGET FILES
        target_files = ["stdout", "stderr", "stdout.submit", "stderr.submit"]
        output = ""
        max_lines = 1000
        for root_dir, subdirs, files in os.walk(run_dir):
            if not root_dir.endswith("execution"):
                continue
            rc_file = os.path.join(root_dir, "rc")
            exitcode = None
            if not os.path.isfile(rc_file):
                continue
            with open(rc_file, "r") as rc_in:
                exitcode = rc_in.readlines()[0].strip()
            if exitcode == "0":
                continue
            output = output + f"[{rc_file}]\n{exitcode}\n" + "-" * 20 + "\n"
            for fname in files:
                if fname not in target_files:
                    continue
                fullname = os.path.join(root_dir, fname)
                lines, is_truncated = tail(fullname, max_lines)
                if is_truncated:
                    output += f"[{fullname}]  (limited to last {max_lines} lines)\n"
                else:
                    output += f"[{fullname}]\n"
                for line in lines:
                    output = output + line
                output += "-" * 20 + "\n"
        return self.success(output)


def tail(filename, max_lines=1000):
    """Return the last n lines of a file.

    :param filename: path to file
    :type filename: str
    :param n: Number of lines to return
    :type n: int
    :return: last n lines of the file and flag if truncated
    :rtype: list[str], bool
    """
    lines = []
    num_lines = 0
    is_truncated = False
    with open(filename, "r") as fh:
        for line in fh:
            lines.append(line)
            num_lines += 1
            if num_lines > max_lines:
                del lines[0]
                is_truncated = True
    return lines, is_truncated


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """JAWS-Site RPC functions"""
    pass


@cli.command()
def server_status():
    """Check Cromwell server status."""
    params = {}
    dispatch = Dispatcher(config.conf)
    result = dispatch.dispatch("server_status", params)
    print(result)


@cli.command()
@click.argument("cromwell_id")
def run_metadata(cromwell_id: str):
    """Get Cromwell metadata for a run."""
    params = {"cromwell_id": cromwell_id}
    dispatch = Dispatcher(config.conf)
    result = dispatch.dispatch("run_metadata", params)
    print(result)


@cli.command()
@click.argument("cromwell_id")
def task_status(cromwell_id: str):
    """Get status of a run."""
    params = {"cromwell_id": cromwell_id}
    dispatch = Dispatcher(config.conf)
    result = dispatch.dispatch("task_status", params)
    print(result)


@cli.command()
@click.argument("cromwell_id")
def task_ids(cromwell_id: str):
    """Retrieve the task IDs for a run."""
    params = {"cromwell_id": cromwell_id}
    dispatch = Dispatcher(config.conf)
    result = dispatch.dispatch("task_ids", params)
    print(result)


@cli.command()
@click.argument("cromwell_id")
def cancel_run(cromwell_id: str):
    """Cancel a run."""
    params = {"cromwell_id": cromwell_id}
    dispatch = Dispatcher(config.conf)
    result = dispatch.dispatch("cancel_run", params)
    print(result)


@cli.command()
@click.argument("cromwell_id")
def run_logs(cromwell_id: str):
    """Get the logs for a run."""
    params = {"cromwell_id": cromwell_id}
    dispatch = Dispatcher(config.conf)
    result = dispatch.dispatch("run_logs", params)
    print(result)


@cli.command()
@click.argument("cromwell_id")
def failure_logs(cromwell_id: str):
    """Get the logs for failed tasks."""
    params = {"cromwell_id": cromwell_id}
    dispatcher = Dispatcher(config.conf)
    result = dispatcher.dispatch("failure_logs", params)
    print(result)
