"""
jaws-site RPC operations for the jaws-workers
"""

from http.client import responses
import logging
import sqlalchemy.exc
from datetime import datetime
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.database import Session
from jaws_site.models import Job_Log
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


def cromwell_submit_job(params):
    """
    Receive a new job submission request (from Cromwell).
    The job is saved in the database and the job_id is returned.
    """
    job_id = None  # TODO
    result = {"job_id": job_id}
    return success(result)

def cromwell_cancel_job(params):
    """
    Receive an abort job command (from Cromwell).
    """
    # TODO
    return success()


def cromwell_is_alive(params):
    """
    Job status check query (from Cromwell).
    The status is retrieved from the db; the worker is not queried.
    """
    # TODO
    is_alive = True
    result = {"is_alive": is_alive}
    return success(result)


def task_complete(params):
    """
    Receive notification of a job's state change (from jaws-worker) and save in db.
    This happens when a job either fails or completes successfully.

    :param params: message from jaws-worker
    :type params: dict
    :return: JSON-RPC2 response
    :rtype: dict
    """
    run_id = params["run_id"]  # JAWS run id
    cromwell_run_id = params["cromwell_run_id"]  # Cromwell's run/workflow UUID
    cromwell_job_id = params["cromwell_job_id"]  # AKA task id
    status_from = params["status_from"]  # always "running"
    status_to = params["status_to"]  # either "complete" or "failed"
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = None
    if "reason" in params:
        reason = params["reason"]
    if reason is not None:
        logger.info(f"Received job status: {cromwell_run_id}:{cromwell_job_id}:{status_from}:{status_to}:{reason}")
    else:
        logger.info(f"Received job status: {cromwell_run_id}:{cromwell_job_id}:{status_from}:{status_to}")

    # DEFINE ROW
    try:
        job_log = Job_Log(
            run_id=run_id,
            cromwell_run_id=cromwell_run_id,
            cromwell_job_id=cromwell_job_id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
        )
    except Exception as error:
        logger.exception(f"Failed to create job_log object for {params}: {error}")
        return failure(500, f"Failed to create job_log object for {params}: {error}")

    # INSERT OR IGNORE
    session = Session()
    try:
        session.add(job_log)
        session.commit()
        logger.debug(f"Run {run_id}, Job {cromwell_job_id} status saved")
    except sqlalchemy.exc.IntegrityError:
        # catch duplicate messages; log event but ignore.  This should never happen.
        session.rollback()
        logger.warning(f"Run {run_id}. Job {cromwell_job_id} status is duplicate; ignored")
    except Exception as error:
        session.rollback()
        session.close()
        logger.exception(f"Run {run_id}. Job {cromwell_job_id} status not saved: {error}")
        return failure(500, f"Failed to insert job {cromwell_job_id} log: {error}")
    session.close()
    return success()


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "update_job_status": {"function": update_job_status, "required_params": [
        "run_id", "cromwell_run_id", "cromwell_job_id", "status_from", "status_to", "timestamp"]}
}
