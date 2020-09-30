"""
jaws-task RPC operations
"""

import logging
from jaws_site import config, api
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def get_task_log(params):
    """Retrieve task log"""
    run_id = params["run_id"]
    try:
        result = api.get_task_log(run_id)
    except api.DatabaseError as error:
        return failure(500, f"Task service db error: {error}")
    return success(result)


def get_task_status(params):
    """
    Retrieve the current status of each task.
    """
    run_id = params["run_id"]
    try:
        result = api.get_task_status(run_id)
    except api.DatabaseError as error:
        return failure(500, f"Task service db error: {error}")
    return success(result)


def update_job_status(params):
    """
    A JTM worker shall post changes in job state, although it is missing the JAWS run id.
    The state change is simply saved in the db; any other actions will be performed by the daemon.
    """
    cromwell_run_id = params["cromwell_run_id"]  # Cromwell's run/workflow UUID
    cromwell_job_id = params["cromwell_job_id"]  # JTM's task_id
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
    "update_job_status": {
        "method": update_job_status,
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
