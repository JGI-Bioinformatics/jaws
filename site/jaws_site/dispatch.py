"""
                20180215.DECKv1b_abrupt4xCO2.ne30_oEC.edison_01_000101_015001_climo.nc
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
import os

from jaws_site.config import jaws_config


logger = logging.getLogger(__package__)


CROMWELL_ENGINE_STATUS_URL = jaws_config.get("CROMWELL", "engine_status_url")
WORKFLOWS_URL = jaws_config.get("CROMWELL", "workflows_url")


def dispatch(method, params):
    """Given a method keyword, call the appropriate function with the provided
    parameters.

    :param method: The method requested.  Abort if unrecognized.
    :type method: str
    :param params: Associated parameters, if any.  Varies by method.
    :type params: dict
    :return: Response in JSON-RPC2 format.
    :rtype: dict
    """
    operations = {
        "server_status": server_status,
        "task_status": task_status,
        "run_metadata": run_metadata,
        "cancel_run": cancel_run,
        "run_logs": run_logs,
        "failure_logs": failure_logs,
    }
    proc = operations.get(method)
    if proc:
        return proc(params)
    return failure(400, "Unknown method")


def success(result):
    """Return a JSON-RPC2 successful result message.

    :param result: The result returned by a successful RPC call.
    :type result: str|int|dict|list
    :return: JSON-RPC2 formatted response
    :rtype: str
    """
    return {"jsonrpc": "2.0", "result": result}


def failure(code, message=None):
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


def do_request(url, method):
    """
    Executes a HTTP method.

    Given the method parameter, it will execute whatever method with the
    provided url. This will usually be a GET or POST.
    """
    try:
        res = method(url)
    except requests.exceptions.RequestException:
        msg = "Cromwell server timeout"
        logger.warning(msg)
        return failure(503, msg)
    return res


def server_status(params):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    response = do_request(CROMWELL_ENGINE_STATUS_URL, requests.get)
    if response.status_code != requests.codes.ok:
        return failure(response.status_code)
    response = {"Cromwell": "UP"}
    return success(response)


def run_metadata(params):
    """Retrieve the metadata of a run.

    :param cromwell_id: Cromwell run ID
    :type params: dict
    :return: The Cromwell metadata for the specified run.
    :rtype: dict
    """

    url = f"{WORKFLOWS_URL}/{params['cromwell_id']}/metadata"

    r = do_request(url, requests.get)
    if r.status_code != requests.codes.ok:
        return failure(r.status_code)
    return success(r.json())


def sort_tasks(metadata, key=None):
    result = []
    for task_name in metadata['calls']:
        task_info = metadata["calls"][task_name]
        job_step_info = task_info[0]
        start = job_step_info.get("start", "?")
        end = job_step_info.get("end", "?")
        execution_status = job_step_info.get("executionStatus")
        result.append((task_name, execution_status, start, end))
    result.sort(key=key)
    return result


def task_status(params):
    """Returns the status for each task processed thus far in a run.

    :param cromwell_id: The Cromwell run ID
    :type cromwell_id: str
    :return: A list of task:status tuples, ordered by start time.
    :rtype: list
    """
    url = f"{WORKFLOWS_URL}/{params['cromwell_id']}/metadata"
    r = do_request(url, requests.get)
    if r.status_code != requests.codes.ok:
        return failure(r.status_code)
    metadata = r.json()

    start_time_key = 2  # our sort key is located in position 2 of the tuple
    result = sort_tasks(metadata, key=lambda x: x[start_time_key])
    return success(result)


def get_tasknames(metadata):
    """
    Helper method that gets all the task names from the cromwell
    job metadata JSON.
    """
    tasks = []
    for task_name in metadata["calls"]:
        task = metadata[task_name]
        if "jobId" in task[-1]:
            job_id = task[-1]
            tasks.append(job_id)
    return tasks


def task_ids(params):
    """Returns a list of all JTM task ids associated with a run thus far.

    :param cromwell_id: The Cromwell run ID
    :type cromwell_id: str
    :return: A list of JTM task IDs.  If a task was rerun,
    only the last is included.
    :rtype: list
    """
    url = f"{WORKFLOWS_URL}/{params['cromwell_id']}/metadata"

    r = do_request(url, requests.get)

    if r.status_code != requests.codes.ok:
        return failure(r.status_code)

    metadata = r.json()
    tasks = get_tasknames(metadata)
    return success(tasks)


def cancel_run(params):
    """Cancel a run.

    :param cromwell_id: The Cromwell run ID
    :type cromwell_id: str
    :return: Either a JSON-RPC2-compliant success or failure message,
    if the run could be cancelled or not.
    :rtype: dict
    """
    url = f"{WORKFLOWS_URL}/{params['cromwell_id']}/abort"
    r = do_request(url, requests.post)
    if r.status_code != requests.codes.ok:
        return failure(r.status_code)
    result = r.json()
    result["status"] = "Cancelling"
    return success(result)


def run_logs(params):
    """
    Retrieve the Cromwell logs for a run.
    """
    url = f'{WORKFLOWS_URL}/{params["cromwell_id"]}/logs'

    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        logger.warning("Cromwell server timeout")
        return failure(500, "Cromwell server timeout")
    if r.status_code != requests.codes.ok:
        return failure(r.status_code)
    return success(r.json())


def rc_exit_code(rc_file):
    """
    Parses the rc file for an exit code
    """
    with open(rc_file, "r") as rc_in:
        exit_code = rc_in.readline().strip()
    return exit_code


def contains_failed_rc(root_dir):
    """
    Checks whether the exit code in an rc file is a non-zero exit code.

    It first checks whether we are in an execution directory and whether
    that execution directory contains an rc file. If it does it will
    check the exit code.
    """
    failed_rc_present = True
    rc = os.path.join(root_dir, "rc")
    if not root_dir.endswith("execution") and not os.path.isfile(rc):
        failed_rc_present = False
    elif rc_exit_code(rc) == "0":
        failed_rc_present = False
    return failed_rc_present


def find_failure_logs(run_dir):
    max_lines = 1000
    separator = "-" * 20 + "\n"
    target_files = ["stdout", "stderr", "stdout.submit", "stderr.submit"]
    output = ""
    for root_dir, subdirs, files in os.walk(run_dir):
        if not contains_failed_rc(root_dir):
            continue
        rc_file = os.path.join(root_dir, "rc")
        output += f"[{rc_file}]\n{rc_exit_code(rc_file)}\n" + separator
        for fname in files:
            if fname not in target_files:
                continue
            full_filepath = os.path.join(root_dir, fname)
            lines, is_truncated = tail(full_filepath, max_lines)
            output += f"[{full_filepath}] "
            if is_truncated:
                output += " (limited to last {max_lines} lines \n"
            else:
                output += "\n"
            for line in lines:
                output += line
            output += " " + separator
    return output


def failure_logs(params):
    """
    Concatenate stdout, stderr files of any failed tasks.
    """
    # GET RUN FOLDER
    url = f'{WORKFLOWS_URL}/{params["cromwell_id"]}/metadata'
    r = do_request(url, requests.get)

    if r.status_code != requests.codes.ok:
        return failure(r.status_code)

    run_dir = r.json()["workflowRoot"]
    output = find_failure_logs(run_dir)
    return success(output)


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
