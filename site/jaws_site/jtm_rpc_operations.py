"""
RPC functions for JTM.
"""

import logging
from datetime import datetime
from jaws_rpc.responses import success, failure
from jaws_site.tasks import TaskLog


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def update_job_status(params, session):
    """JTM shall post changes in job state.  Records contain only cromwell_run and cromwell_job ids, and are
    missing the JAWS run id.  The state change is simply saved in the db."""
    logger.debug(f"{params['cromwell_run_id']} is now {params['status_to']}")
    cromwell_run_id = params["cromwell_run_id"]
    try:
        task_log = TaskLog(session, cromwell_run_id=cromwell_run_id)
        task_log.save_job_log(
            params["cromwell_job_id"],
            params["status_from"],
            params["status_to"],
            datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S"),
            params["reason"],
        )
    except Exception as error:
        logger.error(
            f"Error saving job log for {cromwell_run_id} to {params['status_to']}: {error}"
        )
        return failure(error)
    return success()


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
            "reason",
        ],
    }
}
