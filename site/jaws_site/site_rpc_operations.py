"""
jaws-site has two RPC servers, one on the Central RMQ server and one on the Site RMQ server.
The latter may be on localhost and not accessible from outside.

This package provides the RPC operations for the Site's operations, namely receiving
instructions from Cromwell via jaws-backend.
"""

from jaws_rpc.responses import success, failure
from jaws_site.Task import Task
from jaws_site.Worker import Worker


def worker_begin(params):
    """
    Once a jaws-worker starts, it sends a message to let Site know it's active.  Update db record.
    """
    worker = Worker(worker_id)
    worker.start(params["pid"], params["array_task_id"])


def submit(params):
    """
    Receive a new task submission request.
    The job is saved in the database and the task_id is returned.
    The task_daemon is responsible for actually sending the task to a worker.
    """
    task = Task()
    try:
        task.new(params)
    except Exception as error:
        return failure(500, f"Failed to submit task: {error}")
    return success({"task_id": task.id})


def kill(params):
    """
    Receive an abort task command (from Cromwell).
    """
    try:
        task = Task(params["task_id"])
    except Task.TaskNotFound as error:
        return failure(404, f"Failed to kill task: {error}")
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
    try:
        task = Task(params["task_id"])
    except Task.TaskNotFound as error:
        return failure(404, f"Cannot check Task: {error}")
    return success(task.is_alive)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "worker_begin": {
        "function": worker_begin,
        "required_params": ["worker_id", "array_task_id", "pid"]
    },
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
