"""
jaws-task RPC operations
"""

import logging
from datetime import datetime
from jaws_site import config
from jaws_site.api import Task, TaskNotFoundError, DatabaseError, RpcCommunicatonError, WorkerError
from jaws_rpc.client import RpcClient
from jaws_rpc.responses import success, failure


# logging must be initialized before importing this module
logger = logging.getLogger(__package__)


# used by jaws-worker
def update(params):
    """
    Update a task with new state.

    :param params: cromwell_run_id, cromwell_job_id, status_from, status_to, timetamp, (reason)
    :type params: dict
    :return: JSON-RPC2 response (success/failure); no content
    """
    cromwell_run_id = params["cromwell_run_id"]
    cromwell_job_id = params["cromwell_job_id"]
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = params.get("reason", None)
    logger.info(f"Update {cromwell_run_id}:{cromwell_job_id} status to {status_to}")

    try:
        task = Task(task_id)
        task.update(cromwell_run_id, status_from, status_to, timestamp, reason)
    except TaskNotFoundError:
        return failure(404, "Task not found")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    except ValueError as error:
        return failure(401, f"Invalid task log: {error}")
    return success()


# TODO CHANGE TO RUN_ID
# used by jaws-central
def get_log(params):
    """
    Retrieve all state transitions for each requested task.
 
    :param params: contains "task_ids" list
    :type param: dict
    :return: JSON-RPC2 response; successful result is dict of task to list of logs
    :rtype: dict
    """
    run_id = params["run_id"]

    # get rpc client
    run_service = RunService()

    # get all task_ids associated with the run_id from the jaws-run service
    try:
        task_ids = run_service.get_tasks(run_id)
    except RpcError as error:
        return failure(500, f"RPC failure: {error}; please try again later")
    except RunServiceError as error:
        return failure(500, f"Unable to retrieve task_ids at this time: {error}")

    # get log by task_id
    results = {}
    for task_id in task_ids:
        try:
            task = Task(task_id)
            results[task_id] = task.log()
        except TaskNotFoundError:
            return failure(404, f"Task not found: {task_id}")
        except DatabaseError as error:
            return failure(500, f"Task service db error: {error}")
    return success(results)


# TODO CHANGE TO RUN_ID
# used by jaws-central
def get_status(params):
    """
    Retrieve the current status of each requested task.

    :param params: contains "task_ids"
    :type params: dict
    :return: JSON-RPC2 response; successful result is dict of task_id to status (str)
    :rtype: dict
    """
    run_id = params["run_id"]

    # get rpc client
    run_service = RunService()

    # get all task_ids associated with the run_id from the jaws-run service
    try:
        task_ids = run_service.get_tasks(run_id)
    except RpcError as error:
        return failure(500, f"RPC failure: {error}; please try again later")
    except RunServiceError as error:
        return failure(500, f"Unable to retrieve task_ids at this time: {error}")

    # get status by task_id
    results = {}
    for task_id in task_ids:
        try:
            task = Task(task_id)
            results[task_id] = task.status()
        except TaskNotFoundError:
            return failure(404, f"Task not found: {task_id}")
        except DatabaseError as error:
            return failure(500, f"Task service db error: {error}")
    return success(results)


# used by jaws-run
def any_running(params):
    """
    Given a list of task_ids, determine if any of them are running yet.
    If so, the timestamp of the first task to begin execution shall be returned.
    This is used by the jaws-run service (once the first task is "running",
    the Run will transition from "queued" to "running" state).
    """
    task_ids = params["task_ids"]
    # TODO move this logic to Run API
    first_began = None
    for task_id in task_ids:
        try:
            task = Task(task_id)
            began = task.status()
            if began:
                if not first_began:
                    first_began = began
                elif began < first_began:
                    first_began = began
        except TaskNotFoundError:
            return failure(404, f"Task not found: {task_id}")
        except DatabaseError as error:
            return failure(500, f"Task service db error: {error}")
    result = {
        "running": False,
        "timestamp": None
    }
    if first_began:
        result["running"] = True
        result["timestamp"] = first_began  # TODO string format
    return success(result)


# used by cromwell-backend
def submit(params):
    """
    Submit a task for execution.

    :param params:
    :type params: dict
    :return: JSON-RPC2 response with task_id
    :rtype: dict
    """
    try:
        task = Task(params=params)  # TODO CHANGE TO KVARGS
    except Exception as error:
        return failure(500, f"{error}")
    return success(task.task_id)


# used by cromwell-backend
def check_alive(params):
    """
    Verify the task is running by querying jaws-worker.

    :param params: task_id
    :type params: dict
    :return: JSON-RPC2 response; successful result is True if alive; False otherwise
    :rtype: dict
    """
    task_id = params["task_id"]
    try:
        task = Task(task_id)
        result = task.check_alive()
    except TaskNotFoundError:
        return failure(404, f"Task not found: {task_id}")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    except RpcCommunicationError as error:
        return failure(500, f"RPC communication error: {error}")
    except WorkerError as error:
        return failure(500, f"Worker error: {error}")
    return success(result)


# used by cromwell-backend
def kill(params):
    """
    Abort the execution of a task.

    :param params: task_id
    :type params: dict
    :return: JSON-RPC2 response (success/failure); no content
    :rtype: dict
    """
    task_id = params["task_id"]
    try:
        task = Task(task_id)
        task.cancel()
    except TaskNotFoundError:
        return failure(404, f"Task not found: {task_id}")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    except RpcCommunicationError as error:
        return failure(500, f"Rpc communication error: {error}")
    except WorkerError as error:
        return failure(500, f"Worker error: {error}")
    return success()


# RPC Server dispatch table:
rpc_methods = {
    "update": {
        "method": update,
        "required_params": [
            "cromwell_run_id",
            "cromwell_job_id",
            "status_from",
            "status_to",
            "timestamp",
        ],
    },
    "get_log": {
        "method": get_log,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "get_status": {
        "method": get_status,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "submit": {
        "method": submit,
        "required_params": [
            "script",
            "job-name",
            "cwd",
            "out",
            "err",
            "max-time",
            "memory-gb",
        ],
    },
    "check_alive": {"method": check_alive, "required_params": ["task_id"]},
    "kill": {"method": kill, "required_params": ["task_id"]},
}
