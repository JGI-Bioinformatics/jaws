"""
RPC functions for Central.  Most of these simply wrap the localhost Cromwell REST functions.
"""

import requests
from http.client import responses
import logging
import sqlalchemy.exc
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.database import Session
from jaws_site.models import Run


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


# HELPER FUNCTIONS FOR FORMATTING JSON-RPC2 RESPONSES:


def _success(result={}):
    """Return a JSON-RPC2 successful result message.

    :param result: The result returned by a successful RPC call.
    :type result: str|int|dict|list
    :return: JSON-RPC2 formatted response
    :rtype: str
    """
    return {"jsonrpc": "2.0", "result": result}


def _failure(code, message=None):
    """Return a JSON-RPC2 error message.

    :param code: The error code returned by the RPC call
    :type code: int
    :param message: The error message returned by the RPC call
    :type message: str, optional
    :return: JSON-RPC2 formatted response
    :rtype: dict
    """
    if message is None:
        message = (
            responses["status_code"] if "status_code" in responses else "Unknown error"
        )
    return {"jsonrpc": "2.0", "error": {"code": code, "message": message}}


# RPC OPERATIONS:


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
        return _failure(error.response.status_code, f"Failed to get server status: {error}")
    return _success(status)


def run_metadata(params):
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
        return _success(
            f"Run {run_id} hasn't been submitted to Cromwell, so has no metadata."
        )
    logger.info(f"{user_id} - Run {run_id} - Get metadata")
    try:
        result = cromwell.get_all_metadata(cromwell_run_id)
    except requests.exceptions.HTTPError as error:
        logger.exception(f"Get metadata for {params['run_id']} failed: {error}")
        return _failure(error.response.status_code, f"Failed to get metadata: {error}")
    except Exception as error:
        logger.exception(f"Get metadata for {params['run_id']} failed: {error}")
        return _failure(500, f"Unable to retrieve metadata: {error}")
    return _success(result)


def cancel_run(params):
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
        return _success()

    cromwell_run_id = run.cromwell_run_id
    status = run.status
    run.status = "cancelled"
    try:
        session.commit()
    except Exception as error:
        logger.exception(f"Error updating Run {run_id}: {error}")
        return _failure(500, f"Error updating db record: {error}")
    logger.debug(f"Run {run_id} cancelled")

    # tell Cromwell to cancel the run if it has been submitted to Cromwell already
    if cromwell_run_id and status in ["submitted", "queued", "running"]:
        logger.debug(f"Run {run_id} is {status}: Instructing Cromwell to cancel")
        try:
            cromwell.abort(cromwell_run_id)
        except requests.exceptions.HTTPError as error:
            logger.exception(f"Error aborting cromwell {cromwell_run_id}: {error}")
            return _failure(500, f"Cromwell returned an error: {error}")
        except Exception as error:
            logger.exception(
                f"Unknown error aborting cromwell {cromwell_run_id}: {error}"
            )
            return _failure(500, f"Failed to instruct Cromwell to abort: {error}")
    result = {"cancel": "OK"}
    return _success(result)


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
        return _failure(400, f"Invalid input; {error}")
    session = Session()
    try:
        session.add(run)
        session.commit()
    except Exception as error:
        logger.exception(f"Failed to insert new Run record: {error}")
        return _failure(500, f"Failed to insert new Run record: {error}")
    return _success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "server_status": {"function": server_status},
    "submit": {
        "function": submit,
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
    "run_metadata": {
        "function": run_metadata,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "cancel_run": {
        "function": cancel_run,
        "required_params": ["user_id", "cromwell_run_id"],
    },
}
