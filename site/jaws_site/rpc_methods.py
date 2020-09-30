"""
JAWS Run Service RPC methods.
"""

import logging
from jaws_site.api import Run, Server, RunNotFoundError, DatabaseError, CromwellError
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def server_status(params):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
    try:
        server = Server()
    except CromwellError as error:
        return failure(500, f"Cromwell offline: {error}")
    return success()


def get_status(params):
    """
    Get the current status of a Run.
    """
    try:
        run = Run(params["run_id"])
    except RunNotFoundError as error:
        return failure(404, f"Run not found: {error}")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    return success(run.get_status())


def get_log(params):
    """
    Get the complete logs for a Run.
    """
    try:
        run = Run(params["run_id"])
    except RunNotFoundError as error:
        return failure(404, f"Run not found: {error}")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    return success(run.get_log())


def get_statuses(params):
    """
    Get the current status for several Runs.
    """
    result = {}
    for run_id in params["run_ids"]:
        try:
            run = Run(run_id)
        except RunNotFoundError as error:
            return failure(404, f"Run not found: {error}")
        except DatabaseError as error:
            return failure(500, f"Run db offline: {error}")
        result[run_id] = run.get_status()
    return success(result)


def get_metadata(params):
    """Retrieve the metadata of a run.

    :param run_id: Run ID
    :type params: dict
    :return: The Cromwell metadata for the specified run.
    :rtype: dict
    """
    try:
        run = Run(params["run_id"])
    except RunNotFoundError as error:
        return failure(404, f"Run not found: {error}")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_metadata()
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    return success(result)


def cancel(params):
    """Cancel a run.

    :param cromwell_run_id: The Cromwell run ID
    :type cromwell_run_id: str
    :return: Either a JSON-RPC2-compliant success or failure message,
    :rtype: dict
    """
    try:
        run = Run(params["run_id"])
    except RunNotFoundError as error:
        return failure(404, f"Run not found: {error}")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.cancel()
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    return success(result)


def get_errors(params):
    """Retrieve error messages and stderr for failed Tasks.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: error messages and stderr for failed Tasks
    :rtype: dict
    """
    try:
        run = Run(params["run_id"])
    except RunNotFoundError as error:
        return failure(404, f"Run not found: {error}")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_errors()
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    except IOError as error:
        return failure(500, f"IO error: {error}")
    return success(result)


def submit(params):
    """
    Submit new run
    """
    try:
        run = Run(params["run_id"], params)
    except RunNotFoundError as error:
        return failure(404, f"Run not found: {error}")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    if not run:
        # this doesn't happen, here so flake8 doesn't complain
        return failure(500, "Run object not initialized")
    return success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
rpc_methods = {
    "server_status": {"method": server_status},
    "submit": {
        "method": submit,
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
    "get_status": {"method": get_status, "required_params": ["run_id"]},
    "get_log": {"method": get_log, "required_params": ["run_id"]},
    "get_statuses": {"method": get_statuses, "required_params": ["run_ids"]},
    "get_metadata": {"method": get_metadata, "required_params": ["run_id"],},
    "cancel": {"method": cancel, "required_params": ["run_id"]},
    "get_errors": {"method": get_errors, "required_params": ["run_id"],},
}
