"""
RPC functions for JTM.
"""

from http.client import responses
import logging
import sqlalchemy.exc
from datetime import datetime
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.database import Session
from jaws_site.models import Job_Log


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


def update_job_status(params):
    """JTM shall post changes in job state, although it is missing the JAWS run id.
    The state change is simply saved in the db; any other actions will be performed by the daemon."""
    cromwell_run_id = params["cromwell_run_id"]  # Cromwell's run/workflow UUID
    cromwell_job_id = params["cromwell_job_id"]  # JTM's task_id
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = ""
    if "reason" in params:
        reason = params["reason"]
    logger.info(
        f"Received job status: {cromwell_run_id}:{cromwell_job_id}:{status_from}:{status_to}:{reason}"
    )

    # DEFINE ROW
    try:
        job_log = Job_Log(
            cromwell_run_id=cromwell_run_id,
            cromwell_job_id=cromwell_job_id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
        )
    except Exception as error:
        logger.exception(f"Failed to create job_log object for {params}: {error}")
        return _failure(500, f"Failed to create job_log object for {params}: {error}")

    # INSERT OR IGNORE
    session = Session()
    result = ""
    try:
        session.add(job_log)
        session.commit()
        logger.debug(f"Job {cromwell_job_id} status saved")
    except sqlalchemy.exc.IntegrityError:
        # JTM sometimes sends duplicate messages; ignore
        session.rollback()
        logger.warning(f"Job {cromwell_job_id} status is duplicate; ignored")
        result = "Ignoring duplicate log entry"
    except Exception as error:
        session.rollback()
        session.close()
        logger.exception(f"Job {cromwell_job_id} status not saved: {error}")
        return _failure(500, f"Failed to insert job {cromwell_job_id} log: {error}")
    session.close()
    return _success(result)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "update_job_status": {
        "function": update_job_status,
        "required_params": [
            "cromwell_run_id",
            "cromwell_job_id",
            "status_from",
            "status_to",
            "timestamp",
        ],
    }
}
