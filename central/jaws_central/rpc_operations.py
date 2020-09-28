"""
RPC operations by which Sites update Run and Job status.
"""

import logging
import sqlalchemy.exc
from datetime import datetime
from http.client import responses
from jaws_central.database import Session
from jaws_central.models_sa import Run, Run_Log


logger = logging.getLogger(__package__)


# HELPER FUNCTIONS:
# these are used by RPC operations to format JSON-RPC2 responses


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


def update_run_logs(params):
    """
    Receive a run status update from a Site, insert into run_logs table and update state in runs table.

    :param params: one run log entry
    :type params: dict
    :return: valid JSON-RPC2 response
    :rtype: dict
    """
    site_id = params["site_id"]
    run_id = int(params["run_id"])
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = None
    if "reason" in params:
        reason = params["reason"]
    logger.info(f"Run {run_id} at Site {site_id} now {status_to}")

    # DEFINE ROW OBJ
    try:
        log = Run_Log(
            run_id=run_id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
        )
    except Exception as error:
        err_msg = f"Invalid run_log: {params}; {error}"
        logger.error(err_msg)
        return _failure(400, err_msg)

    # INSERT OR IGNORE INTO RUN-LOGS TABLE
    session = Session()
    try:
        session.add(log)
        session.commit()
    except sqlalchemy.exc.IntegrityError:
        session.rollback()
        logger.warning(f"Run log {run_id}:{status_to} is duplicate; ignored")
        return _success()
    except Exception as error:
        session.rollback()
        err_msg = f"Error inserting run_log: {error}"
        logger.exception(err_msg)
        return _failure(500, err_msg)

    # UPDATE RUNS TABLE
    run = session.query(Run).get(run_id)
    if run.status == status_to:
        # ignore redundant state change (duplicate message)
        session.close()
        return _success()
    run.status = status_to
    run.updated = timestamp
    if status_to == "submitted":
        run.cromwell_run_id = params["cromwell_run_id"]
    elif status_to == "downloading":
        run.download_task_id = params["download_task_id"]
    elif status_to == "succeeded":
        run.result = "succeeded"
    elif status_to == "failed":
        run.result = "failed"
    try:
        session.commit()
    except Exception as error:
        session.rollback()
        logger.exception(f"Failed to update run status: {error}")
        session.close()
        return _failure(500, f"Error inserting log: {error}")
    session.close()
    return _success()


# all RPC operations are defined in this dispatch table
operations = {
    "update_run_logs": {
        "function": update_run_logs,
        "required_parameters": [
            "site_id",
            "run_id",
            "status_from",
            "status_to",
            "timestamp",
        ],
    },
}
