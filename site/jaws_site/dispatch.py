"""
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
import os
import shutil
from jaws_site import config


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
        "delete_run": delete_run,
    }
    proc = operations.get(method)
    logger = logging.getLogger(__package__)
    if proc:
        logger.info(f"RPC: {method} : {params}")
        return proc(params)
    logger.error(f"Unknown RPC method: {method}")
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
        logger = logging.getLogger(__package__)
        logger.warning(msg)
        res = requests.Response()
        res.status_code = 503
    return res


def server_status(params):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    response = do_request(config.conf.get("CROMWELL", "engine_status_url"), requests.get)
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
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{params['cromwell_id']}/metadata"
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
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{params['cromwell_id']}/metadata"
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
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{params['cromwell_id']}/metadata"

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
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{params['cromwell_id']}/abort"
    r = do_request(url, requests.post)
    if r.status_code != 201:
        return failure(r.status_code)
    result = {"cancel": "OK"}
    return success(result)


def run_logs(params):
    """
    Retrieve the Cromwell logs for a run.
    """
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{params['cromwell_id']}/logs"

    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        logger = logging.getLogger(__package__)
        logger.warning("Cromwell server timeout")
        return failure(500, "Cromwell server timeout")
    if r.status_code != requests.codes.ok:
        return failure(r.status_code)
    return success(r.json())


def _find_failure_logs(run_dir):
    max_lines = 1000
    separator = "\n\n" + "-" * 20 + "\n\n"
    target_files = ["rc", "stdout", "stderr", "stdout.submit", "stderr.submit", "script"]
    output = ""
    for root_dir, subdirs, files in os.walk(run_dir):
        if not root_dir.endswith("/execution"):
            continue
        output += f"TASK: {root_dir}\n"

        rc_file = os.path.join(root_dir, "rc")
        stderr_submit_file = os.path.join(root_dir, "stderr.submit")
        if os.path.isfile(rc_file):
            # TASK RAN, CHECK RC
            rc = None
            with open(rc_file, "r") as fh:
                rc = fh.readline().strip()
            if rc == 0:
                continue
        else:
            # TASK DID NOT RUN, CHECK SUBMISSION LOG
            if not os.path.isfile(stderr_submit_file) or os.path.getsize(stderr_submit_file) == 0:
                continue
        # TASK HAS ERROR
        for fname in files:
            if fname not in target_files:
                continue
            full_filepath = os.path.join(root_dir, fname)
            lines, is_truncated = __tail(full_filepath, max_lines)
            if is_truncated:
                output += f"\n[{fname}] showing only last {max_lines} lines\n" + "\n".join(lines)
            else:
                output += f"\n[{fname}]\n" + "\n".join(lines)
        output += separator
    return output


def failure_logs(params):
    """
    Concatenate stdout, stderr files of any failed tasks.
    """
    # GET RUN FOLDER
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{params['cromwell_id']}/metadata"
    r = do_request(url, requests.get)

    if r.status_code != requests.codes.ok:
        return failure(r.status_code)

    run_dir = r.json()["workflowRoot"]
    output = _find_failure_logs(run_dir)
    return success(output)


def __tail(filename, max_lines=1000):
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


def delete_run(params):
    """Delete a run's output.

    :param cromwell_id: The Cromwell run ID
    :type cromwell_id: str
    :return: Either a JSON-RPC2-compliant success or failure message,
    :rtype: dict
    """
    logger = logging.getLogger(__package__)
    if "cromwell_id" not in params:
        return failure(400, "cromwell_id not in params")
    cromwell_id = params["cromwell_id"]
    run_id = params["run_id"]
    url = f"{config.conf.get('CROMWELL', 'workflows_url')}/{cromwell_id}/metadata"
    r = do_request(url, requests.get)
    run_dir = r.json()["workflowRoot"]
    if not run_dir:
        return failure(404, f"workflowRoot not found for {cromwell_id}")
    result = {}
    logger.info(f"Delete output of run_id {run_id} : cromwell_id {cromwell_id} : {run_dir}")
    try:
        # rmtree may take too long, resulting in server timeout error; rename instead
        dest_dir = os.path.join(os.path.dirname(run_dir), f"{os.path.basename(run_dir)}.IGNORE")
        shutil.move(run_dir, dest_dir)
        result["message"] = f"Purged run {run_id} from cache"
    except Exception as error:
        logger.error(f"Failed to purge run {run_id} in {run_dir}: {error}")
        return failure(500, f"Failed to purge run {run_id} in {run_dir}: {error}")
    return success(result)
