"""
RPC operations by which Sites update Run and Job status.
"""

import logging
import sqlalchemy.exc
from datetime import datetime
from http.client import responses
from jaws_central.database import Session
from jaws_central.models_sa import Run, Run_Log, Job_Log


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


def update_run_status(params):
    """
    Receive a run status update and save in RDb.

    :param params: one status update
    :type params: dict
    :return: valid JSON-RPC2 response
    :rtype: dict
    """
    run_id = int(params["run_id"])
    status = params["status"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    logger.info(f"Run {run_id}: now {status}")

    # update Run
    session = Session()
    run = session.query(Run).get(run_id)
    if run.status == status:
        # ignore redundant state change (duplicate message)
        session.close()
        return _success()
    run.status = status
    run.updated = timestamp
    if status == "submitted":
        run.cromwell_run_id = params["cromwell_run_id"]
        logger.info(f"Run {run_id} {status}: {run.cromwell_run_id}")
    elif status == "downloading":
        run.download_task_id = params["download_task_id"]
        logger.info(f"Run {run_id} {status}: {run.download_task_id}")
    elif status == "succeeded":
        run.result = "succeeded"
    elif status == "failed":
        run.result = "failed"
    try:
        session.commit()
    except Exception as error:
        logger.exception(f"Failed to update run status: {error}")
        session.close()
        return _failure(500, f"Error inserting log: {error}")
    session.close()
    return _success()


def update_run_logs(params):
    """
    Receive a set of run status updates and save in RDb.
    Either all updates are saved or none.

    :param params: one or more status updates
    :type params: dict
    :return: valid JSON-RPC2 response
    :rtype: dict
    """
    logs = params["logs"]
    site_id = params["site_id"]
    num_logs = len(logs)
    logger.debug(f"Site {site_id} sent {num_logs} run logs")
    err_msg = None

    for log_entry in logs:
        # INIT VARS
        run_id = int(log_entry[0])
        status_from = log_entry[1]
        status_to = log_entry[2]
        timestamp = datetime.strptime(log_entry[3], "%Y-%m-%d %H:%M:%S")
        reason = log_entry[4]

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
            err_msg = f"Invalid run_log: {log_entry}: {error}"
            logger.error(err_msg)
            continue

        # INSERT OR IGNORE
        session = Session()
        try:
            session.add(log)
            session.commit()
        except sqlalchemy.exc.IntegrityError:
            session.rollback()
            logger.warning(f"Run log {run_id}:{status_to} is duplicate; ignored")
        except Exception as error:
            session.rollback()
            err_msg = f"Error inserting run_log: {error}"
            logger.exception(err_msg)
        session.close()

    if err_msg:
        return _failure(500, err_msg)
    else:
        return _success()


def update_job_logs(params):
    """
    Receive a set of job status updates and save in RDb.
    Either all updates are saved or none.

    :param params: one or more status updates
    :type params: dict
    :return: valid JSON-RPC2 response
    :rtype: dict
    """
    logs = params["logs"]
    site_id = params["site_id"]
    num_logs = len(logs)
    logger.debug(f"Site {site_id} sent {num_logs} job logs")
    err_msg = None

    for log_entry in logs:
        # INIT VARS
        run_id = int(log_entry[0])
        task_name = log_entry[1]
        attempt = int(log_entry[2])
        cromwell_job_id = int(log_entry[3])
        status_from = log_entry[4]
        status_to = log_entry[5]
        timestamp = datetime.strptime(log_entry[6], "%Y-%m-%d %H:%M:%S")
        reason = log_entry[7]

        # DEFINE ROW OBJ
        try:
            log = Job_Log(
                run_id=run_id,
                task_name=task_name,
                attempt=attempt,
                cromwell_job_id=cromwell_job_id,
                status_from=status_from,
                status_to=status_to,
                timestamp=timestamp,
                reason=reason,
            )
        except Exception as error:
            err_msg = f"Invalid job_log: {log_entry} -- {error}"
            logger.error(err_msg)
            continue

        # INSERT OR IGNORE
        session = Session()
        try:
            session.add(log)
            session.commit()
        except sqlalchemy.exc.IntegrityError:
            session.rollback()
            logger.warning(f"Run log {run_id}:{status_to} is duplicate; ignored")
        except Exception as error:
            session.rollback()
            err_msg = f"Error inserting run_log: {error}"
            logger.exception(err_msg)
        session.close()

    if err_msg:
        return _failure(500, err_msg)
    else:
        return _success()


# all RPC operations are defined in this dispatch table
operations = {
    "update_run_status": {
        "function": update_run_status,
        "required_parameters": ["run_id", "status", "timestamp"],
    },
    "update_run_logs": {
        "function": update_run_logs,
        "required_parameters": ["logs", "site_id"],
    },
    "update_job_logs": {
        "function": update_job_logs,
        "required_parameters": ["logs", "site_id"],
    },
}
