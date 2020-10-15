"""
RPC functions for JTM.
"""

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


def update_job_status(params):
    """JTM shall post changes in job state, although it is missing the JAWS run id.
    The state change is simply saved in the db; any other actions will be performed by the daemon."""
    cromwell_run_id = params["cromwell_run_id"]  # Cromwell's run/workflow UUID
    cromwell_job_id = params["cromwell_job_id"]  # JTM's task_id
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = ""
    if "reason" in params and params["reason"] is not None:
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
        return failure(error)

    # INSERT OR IGNORE
    session = Session()
    try:
        session.add(job_log)
        session.commit()
        logger.debug(f"Job {cromwell_job_id} status saved")
    except sqlalchemy.exc.IntegrityError:
        # JTM sometimes sends duplicate messages; ignore
        session.rollback()
        logger.warning(f"Job {cromwell_job_id} status is duplicate; ignored")
    except Exception as error:
        session.rollback()
        session.close()
        return failure(error)
    finally:
        session.close()
    return success(result)


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
