"""
JAWS Run Service RPC methods.
"""

import logging
from jaws_run.api import (
    Run,
    Engine,
    RunNotFoundError,
    RunAlreadyExistsError,
    DatabaseError,
    CromwellError,
    TaskServiceError,
)
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def engine_status(params):
    """
    Return the current status of the Cromwell service.
    """
    logger.info("Check server status")
    try:
        engine = Engine()
    except CromwellError as error:
        return failure(500, f"Cromwell offline: {error}")
    if engine:
        pass  # dummy statement for flake8
    return success()


def get_status(params):
    """
    Get the current status of a Run.

    Required parameters: run_id
    Returns: Current run status (string)
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

    Required parameters: List of run_ids
    Returns: Dict of run_id and it's current status
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

    Required parameters: run_id
    Returns: Table of state transitions.
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

    Required parameters: run_id
    Returns: dict of cromwell_run_id and Cromwell metadata JSON.
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

    Required parameters: run_id
    Returns: None
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

    Required parameters: run_id
    Retruns: dict of failed tasks and their error messages and stderr
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

    Required parameters: run_id, user_id, submission_id, output_endpoint, output_dir
    Returns: None
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

    Required parameters:
    Returns: dict of tasks and logs
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

    Required parameters:
    Returns: dict of tasks and states
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
    "engine_status": {"method": engine_status},
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
