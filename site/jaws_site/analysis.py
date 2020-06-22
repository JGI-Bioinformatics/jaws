"""
RPC functions.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
import os
import collections
from datetime import datetime
from jaws_site import config, wfcopy
from jaws_site.cromwell import Cromwell
from jaws_site.database import Session
from jaws_site.models import Run, Job_Log


MAX_LINES_PER_OUTFILE = 5000


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


# HELPER FUNCTIONS FOR FORMATTING JSON-RPC2 RESPONSES:


def _success(result):
    """Return a JSON-RPC2 successful result message.

    :param result: The result returned by a successful RPC call.
    :type result: str|int|dict|list
    :return: JSON-RPC2 formatted response
    :rtype: str
    """
    return {"jsonrpc": "2.0", "result": result}


def _failure(code, message=None):
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
            responses["status_code"] if "status_code" in responses else "Unknown error"
        )
    return {"jsonrpc": "2.0", "error": {"code": code, "message": message}}


# RPC OPERATIONS:

def server_status(params):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
    try:
        status = cromwell.status()
    except requests.exceptions.HTTPError as error:
        logger.error(f"Failed to get server status: {error}")
        return _failure(error.response.status_code)
    return _success(status)


def run_metadata(params):
    """Retrieve the metadata of a run.

    :param cromwell_id: Cromwell run ID
    :type params: dict
    :return: The Cromwell metadata for the specified run.
    :rtype: dict
    """
    user_id = params["user_id"]
    run_id = params["run_id"]
    cromwell_id = params["cromwell_run_id"]
    logger.info(f"User {user_id}: Get metadata for Run {run_id}")
    if cromwell_id is None:
        return _success(
            f"Run {run_id} hasn't been submitted to Cromwell, so has no metadata."
        )
    logger.info(f"{user_id} - Run {run_id} - Get metadata")
    try:
        metadata = cromwell.get_metadata(params["cromwell_run_id"])
    except requests.exceptions.HTTPError as error:
        logger.error(f"Get metadata for {params['run_id']} failued: {error}")
        return _failure(error.response.status_code)
    except Exception as error:
        logger.exception(f"Get metadata for {params['run_id']} failued: {error}")
        return _failure(f"Unable to retrieve metadata: {error}")
    return _success(metadata.data)


def cancel_run(params):
    """Cancel a run.

    :param cromwell_id: The Cromwell run ID
    :type cromwell_id: str
    :return: Either a JSON-RPC2-compliant success or failure message,
    if the run could be cancelled or not.
    :rtype: dict
    """
    user_id = params["user_id"]
    run_id = params["run_id"]
    cromwell_id = params["cromwell_run_id"]
    logger.info(f"User {user_id}: Cancel run {run_id}")
    try:
        cromwell.abort(cromwell_id)
    except requests.exceptions.HTTPError as error:
        logger.exception(f"Error aborting cromwell {cromwell_id}: {error}")

    # delete run from database because Site only keeps active Run records in db
    try:
        session = Session()
        run = session.query(Run).get(run_id)
        session.delete(run)
        session.commit()
    except Exception as error:
        logger.exception(f"Error deleting Run {run_id}: {error}")

    return _success()


def output(params):
    """
    Get stdout, stderr file output, optionally limited to failed tasks.
    """
    user_id = params["user_id"]
    run_id = params["run_id"]
    logger.info(f"User {user_id}: Get output of Run {run_id}")
    cromwell_id = params["cromwell_run_id"]
    if cromwell_id is None:
        return _success(
            f"Run {run_id} hasn't been submitted to Cromwell, so has no output yet."
        )
    failed_only = False
    if "failed" in params and params["failed"].upper() == "TRUE":
        failed_only = True
    if failed_only:
        logger.info(f"{user_id} - Run {run_id} - Get outfiles, failed tasks")
    else:
        logger.info(f"{user_id} - Run {run_id} - Get outfiles, all tasks")
    try:
        metadata = cromwell.get_metadata(params["cromwell_run_id"])
    except requests.exceptions.HTTPError as error:
        logger.error(f"Failed to get output for Run {run_id}: {error}")
        return _failure(error.response.status_code)
    workflowRoot = metadata.get("workflowRoot")
    logger.debug(f"Find outfiles under {workflowRoot}")
    out_json = _find_outfiles(workflowRoot, failed_only)
    return _success(out_json)


def _find_outfiles(run_dir, failed_only=False):
    """Traverse directory structure to find target files."""
    out_json = collections.defaultdict(dict)
    files = collections.defaultdict(dict)
    rcs = collections.defaultdict(dict)
    target_files = ["stderr", "stderr.submit", "stdout"]

    # get cromwell files for each task, store rc code and target file paths in files dict
    for taskname, filename in wfcopy.get_files(run_dir, delimiter="."):
        basename = os.path.basename(filename)
        if basename == "rc" and failed_only:
            # read file and save return code
            with open(filename) as fh:
                rc = int(fh.read().strip())
            rcs[taskname][basename] = rc
        if basename in target_files:
            files[taskname][basename] = filename

    # for each task, parse target files for msgs, store in out_json.
    for taskname in files:
        if failed_only and not rcs[taskname]["rc"]:
            continue
        for file_type in files[taskname]:
            filename = files[taskname][file_type]
            lines, is_truncated = __tail(filename, MAX_LINES_PER_OUTFILE)
            output = ""
            if is_truncated:
                output += "NOTE: showing only last {MAX_LINES_PER_OUTFILE} lines\n"
                output += "\n".join(lines)
            else:
                output = "\n".join(lines)
            out_json[taskname][file_type] = output

    return dict(out_json)


def __tail(filename, max_lines=1000):
    """Return the last `max_lines` lines of a file.

    :param filename: path to file
    :type filename: str
    :param max_lines: Number of lines to return
    :type max_lines: int
    :return: last lines of the file and flag if truncated
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


def submit(params):
    """Save new run submission in database.  The daemon shall submit to Cromwell after Globus tranfer completes."""
    user_id = params["user_id"]
    run_id = params["run_id"]
    logger.info(f"User {user_id}: Submit new Run {run_id}")
    try:
        run = Run(
            id=run_id,  # pk used by Central
            user_id=user_id,
            email=params["email"],
            transfer_refresh_token=params["transfer_refresh_token"],
            submission_id=params["submission_id"],
            upload_task_id=params["upload_task_id"],
            output_endpoint=params["output_endpoint"],
            output_dir=params["output_dir"],
            status="uploading",
        )
    except Exception as error:
        logger.exception(f"Invalid submit run input; {error}: {params}")
        return _failure(400, "Invalid input; {error}")
    session = Session()
    try:
        session.add(run)
        session.commit()
    except Exception as error:
        logger.exception(f"Failed to insert new Run record: {error}")
        return _failure(500, f"Failed to insert new Run record: {error}")
    return _success()


def update_job_status(params):
    """JTM shall post changes in job status."""
    cromwell_run_id = params["cromwell_run_id"]  # Cromwell's run/workflow UUID
    cromwell_job_id = params["cromwell_job_id"]  # JTM's task_id
    status_from = params["status_from"]
    status_to = params["status_to"]
    if status_from is None:
        status_from = "None"
    if status_to is None:
        status_to = "None"
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = None
    if "reason" in params:
        reason = params["reason"]
    logger.debug(f"Cromwell run:job {cromwell_run_id}:{cromwell_job_id} now {status_to}")

    # CREATE NEW job_log ENTRY
    session = Session()
    try:
        job_log = Job_Log(
            cromwell_run_id=cromwell_run_id,
            cromwell_job_id=cromwell_job_id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
        )
    except Exception as error:
        logger.exception(f"Failed to create job_log object for {params}: {error}")
        return _failure(500, f"Failed to create job_log object for {params}: {error}")
    try:
        session.add(job_log)
    except Exception as error:
        logger.exception(f"Failed to insert job_log: {job_log}: {error}")
        return _failure(500, f"Failed to insert job_log: {job_log}: {error}")
    session.commit()
    return _success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "server_status": {"function": server_status},
    "submit": {
        "function": submit,
        "required_params": [
            "user_id",
            "run_id",
            "email",
            "transfer_refresh_token",
            "submission_id",
            "upload_task_id",
            "output_endpoint",
            "output_dir",
        ],
    },
    "run_metadata": {"function": run_metadata, "required_params": ["cromwell_id"]},
    "cancel_run": {"function": cancel_run, "required_params": ["cromwell_id"]},
    "output": {"function": output, "required_params": ["cromwell_id"]},
    "update_job_status": {"function": update_job_status, "required_params": [
        "cromwell_run_id", "cromwell_job_id", "status_from", "status_to", "timestamp"]}
}
