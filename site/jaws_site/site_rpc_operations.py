"""
jaws-site has two RPC servers, one on the Central RMQ server and one on the Site RMQ server.
This package provides the RPC operations for the Site's operations, namely receiving
instructions from Cromwell.
"""

import logging
import sqlalchemy.exc
from datetime import datetime
from jaws_site import config
from jaws_site.database import Session
from jaws_site.models import Job_Log
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def submit(params):
    """
    Receive a new task submission request (from Cromwell).
    The job is saved in the database and the job_id is returned.
    """
    # insert into db, get task_id
    task_id = None  # TODO
    # submit to worker via rpc
    # TODO
    result = {"job_id": job_id}
    return success(result)

def kill(params):
    """
    Receive an abort task command (from Cromwell).
    """
    # TODO
    return success()


def check_alive(params):
    """
    Job status check query (from Cromwell).
    The status is retrieved from the db; the worker is not queried.
    Unless polling is on, Cromwell only calls this after a restart.
    """
    # TODO
    is_alive = True
    result = {"is_alive": is_alive}
    return success(result)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "submit": {"function": submit, "required_params": []},
    "kill": {"function": kill, "required_params": []},
    "check_alive": {"function": check_alive, "required_params": []},
}
