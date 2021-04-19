"""
RPC functions for Central.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import logging
from sqlalchemy.exc import IntegrityError
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.models import Run
from jaws_rpc.responses import success, failure
from jaws_site.tasks import TaskLog
from jaws_site import errors


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


def server_status(params, session):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
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
    user_id = params["user_id"]
    run_id = params["run_id"]
    cromwell_run_id = params["cromwell_run_id"]
    logger.info(f"User {user_id}: Get metadata for Run {run_id}")
    if cromwell_run_id is None:
        return success(
            f"Run {run_id} hasn't been submitted to Cromwell, so has no metadata."
        )
    logger.info(f"{user_id} - Run {run_id} - Get metadata")
    try:
        result = cromwell.get_all_metadata(cromwell_run_id)
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
    user_id = params["user_id"]
    run_id = params["run_id"]
    logger.info(f"User {user_id}: Cancel run {run_id}")

    # check if run in active runs tables
    try:
        run = session.query(Run).get(run_id)
    except IntegrityError as error:
        logger.exception(f"Run not found: {run_id}: {error}")
    except Exception as error:
        logger.exception(f"Error selecting on runs table: {error}")
    if not run:
        logger.debug(f"Run {run_id} not found")
        return success()

    cromwell_run_id = run.cromwell_run_id
    status = run.status
    run.status = "cancelled"
    try:
        session.commit()
    except Exception as error:
        session.rollback()
        return failure(error)
    logger.debug(f"Run {run_id} cancelled")

    # tell Cromwell to cancel the run if it has been submitted to Cromwell already
    if cromwell_run_id and status in ["submitted", "queued", "running"]:
        logger.debug(f"Run {run_id} is {status}: Instructing Cromwell to cancel")
        try:
            cromwell.abort(cromwell_run_id)
        except Exception as error:
            return failure(error)
    result = {"cancel": "OK"}
    return success(result)


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
        return failure(error)
    return success(errors_report)


def submit(params, session):
    """Save new run submission in database.  The daemon shall submit to Cromwell after Globus tranfer completes."""
    user_id = params["user_id"]
    run_id = params["run_id"]
    logger.info(f"User {user_id}: Submit new Run {run_id}")
    try:
        run = Run(
            id=run_id,  # pk used by Central
            user_id=user_id,
            email=params["email"],
            submission_id=params["submission_id"],
            upload_task_id=params["upload_task_id"],
            output_endpoint=params["output_endpoint"],
            output_dir=params["output_dir"],
            status="uploading",
        )
    except Exception as error:
        return failure(error)
    try:
        session.add(run)
        session.commit()
    except Exception as error:
        session.rollback()
        return failure(error)
    return success()


def get_task_log(params, session):
    """Retrieve task log from database"""
    try:
        tasks = TaskLog(session)
        result = tasks.get_task_log(params["run_id"])
    except Exception as error:
        return failure(error)
    else:
        return success(result)


def get_task_status(params, session):
    """
    Retrieve the current status of each task.
    """
    try:
        tasks = TaskLog(session)
        result = tasks.get_task_status(params["run_id"])
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
    "cancel_run": {
        "function": cancel_run,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "get_task_log": {
        "function": get_task_log,
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
