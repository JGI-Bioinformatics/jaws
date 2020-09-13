"""
jaws-site has two RPC servers, one on the Central RMQ server and one on the Site RMQ server.
The latter may be on localhost and not accessible from outside.

This package provides the RPC operations for the Site's operations, namely receiving
instructions from Cromwell via jaws-backend.
"""

# import logging
from jaws_site import config
from jaws_rpc.responses import success, failure
from jaws_site.Task import Task


# config and logging must be initialized before importing this module
# logger = logging.getLogger(__package__)


def submit(params):
    """
    Receive a new task submission request.
    The job is saved in the database and the task_id is returned.
    The task_daemon is responsible for actually sending the task to a worker.
    """
    try:
        task = Task.new(config.conf, params)
    except Exception as error:
        return failure(500, f"Failed to submit task: {error}")
    return success({"task_id": task.id})


def kill(params):
    """
    Receive an abort task command (from Cromwell).
    """
    task = Task.get(config.conf, params["task_id"])
    try:
        task.kill()
    except Exception as error:
        return failure(500, f"Failed to kill task: {error}")
    return success()


def check_alive(params):
    """
    Check if a Task is running.
    Unless polling is on, Cromwell only calls this after a restart.
    """
    task = Task.get(config.conf, params["task_id"])
    try:
        status = task.status()
    except Exception as error:
        return failure(500, f"Failed to get task status: {error}")
    if status in ("queued", "running"):
        return success(True)
    else:
        return success(False)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "submit": {
        "function": submit,
        "required_params": [
            "script",
            "job_name",
            "cwd",
            "out",
            "err",
            "max_time",
            "memory_gb",
        ],
    },
    "kill": {"function": kill, "required_params": ["task_id"]},
    "check_alive": {"function": check_alive, "required_params": ["task_id"]},
}
