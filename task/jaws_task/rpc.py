"""
jaws-task RPC operations
"""

import logging
from datetime import datetime
from jaws_site import config
from jaws_site.api import Task, TaskNotFoundError, DatabaseError
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
    cromwell_run_id = params["cromwell_run_id"]  # TODO RENAME TO CROMWELL_ID
    cromwell_job_id = params["cromwell_job_id"]  # TODO RENAME TO TASK_ID
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = params.get("reason", None)
    logger.info(f"Job status: {cromwell_run_id}:{cromwell_job_id} is {status_to}")

    try:
        task = Task(task_id)
    except TaskNotFoundError:
        return failure(404, "Task not found")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    try:
        task.update(cromwell_run_id, status_from, status_to, timestamp, reason)
    except ValueError as error:
        return failure(401, f"Invalid task log: {error}")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    return success()


# used by jaws-run
# TODO EITHER BY LIST OF TASK IDS OR BY CROMWELL_RUN_ID
def get_task_log(params):
    """
    Retrieve log of all of a Run's Tasks' state transitions.
 
    :param params: must contain "task_ids" list
    :type param: dict
    :return: JSON-RPC2 response; successful result is table of state transitions
    :rtype: dict
    """
    task_ids = params["task_ids"]
    results = {}
    for task_id in task_ids:
        # TODO
        try:
            log = Log(cromwell_run_id)
        except DatabaseError as error:
            return failure(500, f"Task service db error: {error}")
    return success(results)


# used by jaws-run
# TODO EITHER BY LIST OF TASK IDS OR BY CROMWELL_RUN_ID
def get_task_status(params):
    """
    Retrieve the current status of each Task (so far) of a Run.

    :param params: cromwell_run_id
    :type params: dict
    :return: JSON-RPC2 response; successful result is dict of task_id:status
    :rtype: dict
    """
    cromwell_run_id = params["cromwell_run_id"]
    try:
        log = Log(
            cromwell_run_id
        )  # TODO: NEED TWO CLASSES, ONE FOR RUN_TASKS AND ONE FOR TASK_LOGS
    except DatabaseError as error:
        return failure(500, f"Task service db error: {error}")
    return success(log.get_status())


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
    except TaskNotFoundError:
        return failure(404, f"Task not found: {task_id}")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    worker_id = task.worker_id
    pid = task.pid
    if not worker_id and pid:
        return success(False)

    # request worker status
    rpc_params = config.conf.get_worker_rpc()  # TODO
    rpc_params["queue"] = worker_id
    try:
        worker = RpcClient.client(rpc_params)
    except Exception as error:
        return failure(500, f"{error}")
    response = worker.call("is_alive", {"task_id": task_id, "pid": pid})
    if "error" in response:
        return failure(response["error"]["code"], response["error"]["message"])
    if response["result"]:
        return success(True)
    else:
        return success(False)


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
    except TaskNotFoundError:
        return failure(404, f"Task not found: {task_id}")
    except DatabaseError as error:
        return failure(500, f"Task db error: {error}")
    return success(task.cancel())


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
    "get_task_log": {
        "method": get_task_log,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "get_task_status": {
        "method": get_task_status,
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
