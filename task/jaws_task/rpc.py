"""
jaws-task RPC operations
"""

import logging
from jaws_site import config, api
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def submit(params):
    """
    Submit a task for execution.
    """
    pass


def status(params):
    """
    Check the status of a task.
    """
    pass

def cancel(params):
    """
    Abort the execution of a task.
    """
    pass


def get_task_log(params):
    """
    Retrieve log of all tasks' state transitions.
    
    Required parameters: run_id
    Returns: Table of task state transitions
    """
    run_id = params["run_id"]
    try:
        result = api.get_task_log(run_id)
    except api.DatabaseError as error:
        return failure(500, f"Task service db error: {error}")
    return success(result)


def get_task_status(params):
    """
    Retrieve the current status of each Task of a Run.
    Tasks which haven't been scheduled yet are not shown.

    Required parameters: run_id
    Returns: Current status for each active and completed task.
    """
    run_id = params["run_id"]
    try:
        result = api.get_task_status(run_id)
    except api.DatabaseError as error:
        return failure(500, f"Task service db error: {error}")
    return success(result)


def update_task_status(params):
    """
    Receive a notification from the backend of a task's change in state.

    Required parameters: cromwell_run_id, cromwell_job_id, status_from, status_to, timetamp
    Optional parameters: reason (i.e. error message)
    """
    # A JTM worker shall post changes in job state, but is missing run_id, task_name, and attempt fields.
    cromwell_run_id = params["cromwell_run_id"]
    cromwell_job_id = params["cromwell_job_id"]
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = params.get("reason", None)
    logger.info(f"Job status: {cromwell_run_id}:{cromwell_job_id} is {status_to}")

    try:
        api.update_job_status(cromwell_run_id, cromwell_job_id, status_from, status_to, timestamp, reason)
    except api.DatabaseError as error:
        return failure(500, f"Task service db error: {error}")
    return success()


# RPC Server dispatch table:
rpc_methods = {
    "update_task_status": {
        "method": update_task_status,
        "required_params": [
            "cromwell_run_id",
            "cromwell_job_id",
            "status_from",
            "status_to",
            "timestamp",
        ],
    },
    "get_task_log": {
        "method": get_task_log,
        "required_params": ["user_id", "run_id"],
    },
    "get_task_status": {
        "method": get_task_status,
        "required_params": ["user_id", "run_id"],
    },
}
