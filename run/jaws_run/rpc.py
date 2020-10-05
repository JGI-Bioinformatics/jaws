"""
JAWS Run Service RPC methods.
"""

import logging
from jaws_run.api import (
    Run,
    RunNotFoundError,
    RunAlreadyExistsError,
    DatabaseError,
    CromwellError,
    TaskServiceError,
)
from jaws_run.config import conf
from jaws_run.cromwell import Cromwell
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)
cromwell = Cromwell(conf.get("CROMWELL", "url"))


def get_cromwell_status(params):
    """
    Return the current status of the Cromwell service.

    :param params: None
    :type params: dict
    :return: JSON-RPC2 successful response with no content
    :rtype: dict
    """
    logger.debug("Check Cromwell status")
    try:
        _ = cromwell.status()
    except Exception:
        return failure(500, "Cromwell offline")
    return success()


def get_status(params):
    """
    Get the current status of a Run.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response, successful result is current status
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_status()
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)


def get_statuses(params):
    """
    Get the current status for a list of Runs.

    :param params: List of run_ids
    :type params: dict
    :return: JSON-RPC2 response with result of run_id:status
    :rtype: dict
    """
    results = {}
    for run_id in params["run_ids"]:
        try:
            run = Run(run_id)
        except RunNotFoundError:
            return failure(404, f"Run {run_id} not found")
        except DatabaseError as error:
            return failure(500, f"Run db offline: {error}")
        try:
            result = run.get_status()
        except DatabaseError as error:
            return failure(500, f"Run db offline: {error}")
        except Exception as error:
            return failure(500, f"{error}")
        results.append(result)
    return success(results)


def get_log(params):
    """
    Get a log of a state transitions of a Run.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with table of Run state transitions.
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_log()
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)


def get_metadata(params):
    """
    Retrieve the Cromwell metadata of a run and all it's subworkflows.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with dict of cromwell_run_id:Cromwell metadata JSON.
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_metadata()
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)


def cancel(params):
    """
    Cancel a run.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with no content
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        run.cancel()
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    return success()


def get_errors(params):
    """
    Retrieve error messages and stderr for failed Tasks.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with dict of failed tasks: err msg and stderr
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
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
    Submit new run for execution.

    :param params: run_id, user_id, submission_id, output_endpoint, output_dir
    :type params: dict
    :return: JSON-RPC2 response to no content
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id, params)
    except RunAlreadyExistsError:
        return failure(404, f"Run {run_id} already exists")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    except ValueError as error:
        return failure(401, f"Invalid Run parameters: {error}")
    if run:
        pass  # so flake8 doesn't complain
    return success()


def get_task_log(params):
    """
    Get log of Task state transitions for a Run.
    The tasks are retrieved from Cromwell metadata and state transitions retrieved from jaws-task.

    :param params:
    :type params: dict
    :return: JSON-RPC2 response with dict of tasks and logs
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_task_log()
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    except TaskServiceError as error:
        return failure(500, f"Task service error: {error}")
    return success(result)


def get_task_status(params):
    """
    Get state of each task in a Run.
    The tasks are retrieved from Cromwell metadata and current state retrieved from jaws-task.

    :param params:
    :type params: dict
    :return: JSON-RPC2 response with result dict of task:status
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    try:
        result = run.get_task_status()
    except CromwellError as error:
        return failure(500, f"Cromwell error: {error}")
    except TaskServiceError as error:
        return failure(500, f"Task service error: {error}")
    return success(result)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
rpc_methods = {
    "get_cromwell_status": {"method": get_cromwell_status},
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
    "get_statuses": {"method": get_statuses, "required_params": ["run_ids"]},
    "get_log": {"method": get_log, "required_params": ["run_id"]},
    "get_metadata": {"method": get_metadata, "required_params": ["run_id"]},
    "cancel": {"method": cancel, "required_params": ["run_id"]},
    "get_errors": {"method": get_errors, "required_params": ["run_id"]},
    "get_task_log": {"method": get_task_log, "required_params": ["run_id"]},
    "get_task_status": {"method": get_task_status, "required_params": ["run_id"]},
}
