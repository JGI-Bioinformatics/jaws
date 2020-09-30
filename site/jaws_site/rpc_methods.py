"""
RPC functions for Central.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
import logging
import sqlalchemy.exc

# from sqlalchemy.exc import SQLAlchemyError
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.db import Session, Run
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


def server_status(params):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
    try:
        status = cromwell.status()
    except requests.exceptions.HTTPError as error:
        logger.exception(f"Failed to get server status: {error}")
        return failure(
            error.response.status_code, f"Failed to get server status: {error}"
        )
    return success(status)


def get_status(params):
    """
    Get the current status of a Run.
    """
    try:
        run = Run(params["run_id"])
    except Exception as error:
        return failure(f"Unable to retrieve status at this time: {error}")
    return success(run.get_status())


def get_log(params):
    """
    Get the complete logs for a Run.
    """
    try:
        run = Run(params["run_id"])
    except Exception as error:
        return failure(f"Unable to retrieve status at this time: {error}")
    return success(run.get_log())


def get_statuses(params):
    """
    Get the current status for several Runs.
    """
    result = {}
    for run_id in params["run_ids"]:
        try:
            run = Run(run_id)
        except Exception as error:
            return failure(f"Unable to retrieve status at this time: {error}")
        result[run_id] = run.get_status()
    return success(result)


def get_metadata(params):
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
    except requests.exceptions.HTTPError as error:
        logger.exception(f"Get metadata for {params['run_id']} failed: {error}")
        return failure(error.response.status_code, f"Failed to get metadata: {error}")
    except Exception as error:
        logger.exception(f"Get metadata for {params['run_id']} failed: {error}")
        return failure(500, f"Unable to retrieve metadata: {error}")
    return success(result)


def cancel(params):
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
        session = Session()
        run = session.query(Run).get(run_id)
    except sqlalchemy.exc.IntegrityError as error:
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
        logger.exception(f"Error updating Run {run_id}: {error}")
        return failure(500, f"Error updating db record: {error}")
    logger.debug(f"Run {run_id} cancelled")

    # tell Cromwell to cancel the run if it has been submitted to Cromwell already
    if cromwell_run_id and status in ["submitted", "queued", "running"]:
        logger.debug(f"Run {run_id} is {status}: Instructing Cromwell to cancel")
        try:
            cromwell.abort(cromwell_run_id)
        except requests.exceptions.HTTPError as error:
            logger.exception(f"Error aborting cromwell {cromwell_run_id}: {error}")
            return failure(500, f"Cromwell returned an error: {error}")
        except Exception as error:
            logger.exception(
                f"Unknown error aborting cromwell {cromwell_run_id}: {error}"
            )
            return failure(500, f"Failed to instruct Cromwell to abort: {error}")
    result = {"cancel": "OK"}
    return success(result)


def get_errors(params):
    """Retrieve error messages and stderr for failed Tasks.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: error messages and stderr for failed Tasks
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
        metadata = cromwell.get_metadata(cromwell_run_id)
        result = metadata.errors()
    except requests.exceptions.HTTPError as error:
        logger.exception(f"Get errors for {params['run_id']} failed: {error}")
        return failure(error.response.status_code, f"Failed to get errors: {error}")
    except Exception as error:
        logger.exception(f"Get errors for {params['run_id']} failed: {error}")
        return failure(500, f"Unable to retrieve errors: {error}")
    return success(result)


def submit(params):
    """Save new run submission in database.  The daemon shall submit to Cromwell after Globus tranfer completes."""
    user_id = params["user_id"]
    run_id = params["run_id"]
    logger.info(f"User {user_id}: Submit new Run {run_id}")
    try:
        run = Run(
            id=run_id,  # pk used by Central
            user_id=user_id,
            email=params["email"],
            transfer_refresh_token=params["transfer_refresh_token"],
            submission_id=params["submission_id"],
            upload_task_id=params["upload_task_id"],
            output_endpoint=params["output_endpoint"],
            output_dir=params["output_dir"],
            status="uploading",
        )
    except Exception as error:
        logger.exception(f"Invalid submit run input; {error}: {params}")
        return failure(400, f"Invalid input; {error}")
    session = Session()
    try:
        session.add(run)
        session.commit()
    except Exception as error:
        session.rollback()
        logger.exception(f"Failed to insert new Run record: {error}")
        return failure(500, f"Failed to insert new Run record: {error}")
    return success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
rpc_methods = {
    "server_status": {"method": server_status},
    "submit": {
        "method": submit,
        "required_params": [
            "user_id",
            "run_id",
            "email",
            "transfer_refresh_token",
            "submission_id",
            "upload_task_id",
            "output_endpoint",
            "output_dir",
        ],
    },
    "get_status": {
        "method": get_status,
        "required_params": ["run_id"]
    },
    "get_log": {
        "method": get_log,
        "required_params": ["run_id"]
    },
    "get_statuses": {
        "method": get_statuses,
        "required_params": ["run_ids"]
    },
    "get_metadata": {
        "method": get_metadata,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "cancel": {
        "method": cancel,
        "required_params": ["user_id", "cromwell_run_id"]
    },
    "get_errors": {
        "method": get_errors,
        "required_params": ["user_id", "cromwell_run_id"],
    },
}
