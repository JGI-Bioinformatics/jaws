"""
jaws-worker RPC operations.
"""

import psutil
import os.path
import os.environ
import subprocess, signal, os, threading, errno
from contextlib import contextmanager
from jaws_worker import config
from jaws_worker.database import Session
from jaws_worker.models import Jobs, Job_Log
from jaws_rpc.responses import success, failure

# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)

# array_task_id = os.environ["ARRAY_TASK_ID"]


def quit(params):
    """
    Quit all processes to release hardware reservation.
    """
    quit = threading.Thread(target=__quit)
    quit.start()
    return success()


def __quit():
    """Wait a moment in order to allow RPC response to be returned and terminate all processes."""
    sleep(2)
    pid = os.getpid()
    os.kill(pid, signal.STOP)


def run(params):
    """
    Execute the provided shell script and return the process id.

    :params param: message from jaws-site
    :type params: dict
    :return: JSON-RPC2 response
    :rtype: dict
    """
    cmd = params["cmd"].split()
    try:
        proc = subprocess.Popen(cmd)
    except FileNotFound as error:
        return failure(404, f"FileNotFound: {error}")
    result = {"pid": proc.pid}
    return success(result)


def status(params):
    """
    Check if jaws-worker is alive and able to accept RMQ messages.
    If optional pid provided, check if the Task is running.
    If it's not running, check if 'rc' file exists.
    If so, read rc, else write rc of 1 so that Cromwell can pick it up.
    """
    if "pid" not in params:
        # yes, worker is alive and responsive
        return success()

    # check on a Task
    task_pid = params["pid"]
    is_alive = None
    rc = None
    try:
        task = psutil.Process(task_pid)
    except psutil.NoSuchProcess:
        # process is not running
        is_alive = False

        # check if rc file exists
        if os.path.isfile(rc_file):
            # File exists, process exited cleanly.
            # Read rc value from file
            with open(rc_file, "r") as fh:
                rc = fh.read()
        else:
            # This should not happen unless killed -9.
            # Write rc file with failure code for Cromwell.
            rc = "1"
            with open(rc_file, "w") as fh:
                fh.write(rc)
    else:
        is_alive = True

    result = {"is_alive": is_alive, "rc": rc}
    return success(result)


def check_alive(params):
    """
    Simplified call for checking on a task; rpc_server verifies pid is provided.
    :return: True if alive, False otherwise
    :rtype: bool
    """
    response = status(params)
    if "result" in response:
        return success(response["result"]["is_alive"])
    else:  # error
        return response


def kill(params):
    """
    Recursively kill process and write return code to file.
    Cromwell reads the rc_file to determine status.
    """
    task_pid = params["pid"]
    call_root_dir = params["call_root"]
    rc_file = os.path.join(call_root, "execution", "rc")
    rc = 9  # any nonzero value will do
    if "rc" in params:
        rc = params["rc"]

    # check if running
    try:
        task = psutil.Process(task_pid)
    except psutil.NoSuchProcess:
        # not running; write rc file if it doesn't exist
        if not os.path.isfile(rc_file):
            with open(rc_file, "w") as fh:
                fh.write(str(rc))
    else:
        # recursively kill
        for child in task.children(recursive=True):
            child.kill()
        task.kill()

        # write rc file, if necessary
        if not os.path.isfile(rc_file):
            with open(rc_file, "w") as fh:
                fh.write(str(rc))

    return success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "status": {"function": status},
    "quit": {"function": quit, "required_params": ["task_id"]},
    "run": {"function": run, "required_params": ["cmd"],},
    "check_alive": {"function": check_alive, "required_params": ["pid"]},
    "kill": {"function": kill, "required_params": ["pid"]},
}
