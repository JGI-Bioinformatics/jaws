"""
RPC operations by which Sites update Run and Job status.
"""

import logging
import sqlalchemy.exc
from datetime import datetime
from jaws_central.models import Run, Run_Log, User
from jaws_central.transfers import Transfer
from jaws_rpc.responses import success, failure
from sqlalchemy.exc import SQLAlchemyError

logger = logging.getLogger(__package__)


def update_run_logs(params, session):
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
        return failure(error)

    # INSERT OR IGNORE INTO RUN-LOGS TABLE
    try:
        session.add(log)
        session.commit()
    except sqlalchemy.exc.IntegrityError:
        session.rollback()
        logger.warning(f"Run log {run_id}:{status_to} is duplicate; ignored")
        return success()
    except Exception as error:
        session.rollback()
        return failure(error)

    # UPDATE RUNS TABLE
    run = session.query(Run).get(run_id)
    if run.status == status_to:
        # ignore redundant state change (duplicate message)
        return success()
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
        return failure(error)

    return success()


def _get_email_address(session, user_id):
    """
    Get the email address associated with a given user id.  Raise if not found.
    :param session: sqlalchemy currently db session
    :type session:  sqlalchemy.session
    :param user_id: user id
    :type user_id: str
    :return: email address
    :rtype: str
    """
    assert user_id
    try:
        user = session.query(User).filter(User.id == user_id).one_or_none()
    except SQLAlchemyError as e:
        logger.error(f"Db error: {e}")
    if user is None:
        logger.error(f"Could not found user {user_id}")
    return user.email


def submit_transfer(params, session):
    try:
        transfer = Transfer.from_params(session, params)
    except Exception as error:
        return failure(str(error))
    else:
        return success({"transfer_id": transfer.data.id})


def transfer_status(params, session):
    """
    Request transfer status from specified site (via RPC) and return the results.
    """
    transfer_id = params["transfer_id"]
    try:
        transfer = Transfer.from_id(transfer_id)
    except Exception as error:
        return failure(str(error))
    else:
        return success({"status": transfer.status})


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
    "submit_transfer": {
        "function": submit_transfer,
        "required_parameters": [
            "src_site_id",
            "dest_site_id",
            "manifest_json",
        ],
    },
    "transfer_status": {
        "function": transfer_status,
        "required_parameters": [
            "transfer_id",
        ],
    },
}
