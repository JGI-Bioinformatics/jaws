"""
RPC operations by which Sites update Run and Job status.
"""

import logging
import sqlalchemy.exc
from datetime import datetime
from jaws_central.models_sa import Run, Run_Log
from jaws_rpc.responses import success, failure
from jaws_central.xfer_queue import XferQueue


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
        run.download_id = params["download_id"]
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


def transfer_status(params, session):
    """Query XferQueue for status of a transfer"""
    try:
        xq = XferQueue(session)
        status = xq.transfer_status(params["xfer_id"])
    except Exception as error:
        return failure(error)
    else:
        return success(status)


def submit_transfer(params, session):
    """Submit transfer task to XferQueue and return xfer_id"""
    try:
        xq = XferQueue(session)
        xfer_id = xq.submit_transfer(**params)
    except Exception as error:
        return failure(error)
    else:
        return success(xfer_id)


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
    "transfer_status": {
        "function": transfer_status,
        "required_parameters": [
            "xfer_id",
        ],
    },
    "submit_transfer": {
        "function": submit_transfer,
        "required_parameters": [
            "src_endpoint",
            "dest_endpoint",
            "manifest",
            "label",
            "user",
        ],
    },
}
