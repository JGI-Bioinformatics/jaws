"""
RPC operations used by jaws-sites.
"""

import logging
from datetime import datetime
from jaws_central.runs import Run
from jaws_rpc.responses import success, failure


logger = logging.getLogger(__package__)


def update_run_logs(params, session):
    """
    Receive a run status update from a Site.

    :param params: one run log entry
    :type params: dict
    :return: valid JSON-RPC2 response
    :rtype: dict
    """
    site_id = params["site_id"]
    run_id = int(params["run_id"])
    status_to = params["status_to"]
    params["timestamp"] = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    logger.info(f"Run {run_id} at Site {site_id} now {status_to}")
    try:
        run = Run(session, run_id)
        run.update_status(**params)
    except Exception as error:
        return failure(error)
    return success()


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
            "reason",
        ],
    },
}
