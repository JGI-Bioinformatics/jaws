"""
jaws-worker RPC operations
"""

import logging
import sqlalchemy.exc
from datetime import datetime
from jaws_worker import config
from jaws_worker.database import Session
from jaws_worker.models import Jobs, Job_Log
from jaws_rpc.response import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def health(params):
    """
    Worker health check.

    :param params: message from jaws-site
    :type params: dict
    :return: JSON-RPC2 response
    :rtype: dict
    """
    return success()

def job_status(params):
    """
    Check a job's status.

    :param params: message from jaws-site
    :type params: dict
    :return: JSON-RPC2 response
    :rtype: dict
    """
    job_id = params["job_id"]
    # TODO
    result = True
    return success(result)


def run_job(params):
    """
    Run a job.

    :params param: message from jaws-site
    :type params: dict
    :return: JSON-RPC2 response
    :rtype: dict
    """
    # TODO
    return success()


def cancel_job(params):
    """
    Abort a job.
    """
    # TODO
    return success()


def quit(params):
    """
    Worker should shut down.
    """
    return success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "health": {"function": health, "required_params": []},
    "job_status": {"function": job_status, "required_params": ["job_id"]},
    "run_job": {
        "function": run_job,
        "required_params": ["run_id", "cromwell_run_id", "job_id", "task_name", "attempt", "script", "num_cores", "ram_gb", "minutes"]
    },
    "cancel_job": {"function": cancel_job, "required_params": ["job_id"]},
    "quit": {"function": quit, "required_params": ["job_id"]},
}
