"""
RPC functions for Central.
"""

import logging
from jaws_rpc.responses import success, failure
from jaws_site.tasks import TaskLog
from jaws_site import errors
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.runs import Run, RunNotFound


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def server_status(params, session):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    try:
        status = cromwell.status()
    except Exception as error:
        return failure(error)
    return success(status)


def run_metadata(params, session):
    """Retrieve the metadata of a run.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: The Cromwell metadata for the specified run.
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Metadata Run {params['run_id']}")
    try:
        run = Run(session, **params)
        result = run.metadata()
    except Exception as error:
        return failure(error)
    return success(result)


def run_outputs(params, session):
    """Retrieve the outputs-json of a run.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: The Cromwell outputs for the specified run.
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Outputs Run {params['run_id']}")
    try:
        run = Run(session, **params)
        result = run.outputs()
    except Exception as error:
        return failure(error)
    return success(result)


def cancel_run(params, session):
    """Cancel a run.

    :param cromwell_run_id: The Cromwell run ID
    :type cromwell_run_id: str
    :return: Either a JSON-RPC2-compliant success or failure message,
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Cancel Run {params['run_id']}")
    try:
        run = Run(session, **params)
        run.cancel()
    except RunNotFound as error:
        return success(f"Cancelled with Cromwell warning: {error}")
    except Exception as error:
        return failure(error)
    return success(True)


def get_errors(params, session):
    """Retrieve error report which includes errors from both Cromwell metadata and the TaskLog.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: errors report
    :rtype: dict
    """
    user_id = params["user_id"]
    run_id = params["run_id"]
    cromwell_run_id = params["cromwell_run_id"]
    logger.info(f"User {user_id}: Get errors for Run {run_id}")
    if cromwell_run_id is None:
        return success(f"Run {run_id} hasn't been submitted to Cromwell.")
    logger.info(f"{user_id} - Run {run_id} - Get errors")
    try:
        errors_report = errors.get_errors(session, cromwell_run_id)
    except Exception as error:
        logger.error(f"Error getting errors report for {run_id}: {error}")
        return failure(error)
    return success(errors_report)


def submit(params, session):
    """Save new run submission in database.  The daemon shall submit to Cromwell after Globus tranfer completes."""
    logger.info(f"User {params['user_id']}: Submit Run {params['run_id']}")
    try:
        run = Run(session, **params)
    except Exception as error:
        return failure(error)
    else:
        return success(run.status)


def get_task_log(params, session):
    """Retrieve task log from database"""
    logger.info(f"User {params['user_id']}: Task-log Run {params['run_id']}")
    try:
        task_log = TaskLog(session, run_id=params["run_id"])
        result = task_log.task_log()
    except Exception as error:
        return failure(error)
    else:
        return success(result)


def get_task_summary(params, session):
    """Retrieve task summary from database"""
    logger.info(f"User {params['user_id']}: Task-summary Run {params['run_id']}")
    try:
        task_log = TaskLog(session, run_id=params["run_id"])
        result = task_log.task_summary()
    except Exception as error:
        return failure(error)
    else:
        return success(result)


def get_task_status(params, session):
    """
    Retrieve the current status of each task.
    """
    logger.info(f"User {params['user_id']}: Task-status Run {params['run_id']}")
    try:
        task_log = TaskLog(session, run_id=params["run_id"])
        result = task_log.task_status()
    except Exception as error:
        return failure(error)
    else:
        return success(result)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "server_status": {"function": server_status},
    "submit": {
        "function": submit,
        "required_params": [
            "user_id",
            "run_id",
            "email",
            "submission_id",
            "upload_task_id",
            "output_endpoint",
            "output_dir",
        ],
    },
    "run_metadata": {
        "function": run_metadata,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_outputs": {
        "function": run_outputs,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "cancel_run": {
        "function": cancel_run,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "get_task_log": {
        "function": get_task_log,
        "required_params": ["user_id", "run_id"],
    },
    "get_task_summary": {
        "function": get_task_summary,
        "required_params": ["user_id", "run_id"],
    },
    "get_task_status": {
        "function": get_task_status,
        "required_params": ["user_id", "run_id"],
    },
    "get_errors": {
        "function": get_errors,
        "required_params": ["user_id", "cromwell_run_id"],
    },
}
