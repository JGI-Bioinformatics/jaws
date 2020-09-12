"""
jaws-worker RPC operations.
"""

import subprocess
import psutil
import os
import os.path
import signal
from jaws_worker import config
from jaws_worker.database import Session
from jaws_worker.models import Jobs, Job_Log
from jaws_rpc.response import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def health(params):
    """
    Verify jaws-worker is running and can accept RPC calls.

    :param params: message from jaws-site
    :type params: dict
    :return: JSON-RPC2 response
    :rtype: dict
    """
    return success()


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
    # init vars
    cmd = params["cmd"].split()
    call_root_dir = params["call_root"]
    rc_file = os.path.join(call_root, "execution", "rc")
    max_time = params["time"].split(":")  # format is "hh:mm:ss"
    for n in range(len(max_time), 3):
        max_time.insert(0, 0)  # add any missing fields
    max_sec = max_time[0] * 3600 + max_time[1] * 60 + max_time[2]

    # Run the Task in the background.
    # The worker doesn't care about the Task's exit status or output.
    # jaws-site will know when the Task is complete and won't send the
    # worker another Task before then (or will instruct worker to quit).
    try:
        proc = subprocess.Popen(cmd)
    except FileNotFound as error:
        return failure(404, f"FileNotFound: {error}")

    # return pid
    result = {"pid": proc.pid}
    return success(result)


def check_alive(params):
    """
    Check if Task is running.
    If it's not running, check if 'rc' file exists.
    If so, read rc, else write rc of 1 so that Cromwell can pick it up.
    """
    task_pid = params["task_id"]
    is_alive = True
    rc = None
    try:
        task = psutil.Process(task_pid)
    except psutil.NoSuchProcess:
        # not running
        is_alive = False
        # check if rc file exists
        if os.path.isfile(rc_file):
            with open(rc_file, "r") as fh:
                rc = fh.read()
        else:
            rc = "1"
            with open(rc_file, "w") as fh:
                fh.write(rc)
    result = {"is_alive": is_alive, "rc": rc}
    return success(result)


def kill(params):
    """
    Recursively kill process and write return code to file.
    Cromwell reads the rc_file to determine status.
    """
    task_pid = params["task_id"]
    call_root_dir = params["call_root"]
    rc_file = os.path.join(call_root, "execution", "rc")
    # rc_file=params["rc_file"]
    rc = 9
    if "rc" in params:
        rc = params["rc"]

    try:
        task = psutil.Process(task_pid)
    except psutil.NoSuchProcess:
        # not running; write rc file if it doesn't exist
        if not os.path.isfile(rc_file):
            with open(rc_file, "w") as fh:
                fh.write(str(rc))
        return success()

    # recursively kill
    for child in task.children(recursive=True):
        child.kill()
    task.kill()

    # write rc file
    with open(rc_file, "w") as fh:
        fh.write(str(rc))
    return success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "health": {"function": health, "required_params": []},
    "quit": {"function": quit, "required_params": ["task_id"]},
    "run": {
        "function": run_task,
        "required_params": [
            "run_id",
            "cromwell_run_id",
            "task_id",
            "task_name",
            "attempt",
            "script",
            "num_cores",
            "ram_gb",
            "minutes",
        ],
    },
    "check_alive": {"function": check_alive, "required_params": ["task_id"]},
    "kill": {"function": kill, "required_params": ["task_id"]},
}
