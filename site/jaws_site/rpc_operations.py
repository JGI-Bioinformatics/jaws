import logging
from jaws_rpc.responses import success, failure
from jaws_site import config
from jaws_site import queue_wait as slurm_queue_wait
from jaws_site.cromwell import Cromwell
from jaws_site.runs import Run
from jaws_site.tasks import TaskLog
from jaws_site.transfers import Transfer


DEFAULT_TZ = "America/Los_Angeles"


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


def server_status(params, session):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    try:
        status = cromwell.status()
    except Exception as error:
        return failure(error)
    return success(status)


def queue_wait(params, session):
    """Return the current queue wait times of the possible
    condor pools (sm, md, lg, xlg).

    :return: returns the estimated queue wait times for each condor pool.
    :rtype: dict
    """

    logger.info("Check estimated queue wait times for condor slurm pools")
    try:
        result = slurm_queue_wait.check_queue_wait(logger)
    except Exception as error:
        return failure(error)
    return success(result)


def cancel_run(params, session):
    """Mark a Run to be cancelled.  It will be cancelled by the Run daemon later (asynchronously).

    :param run_id: The JAWS run ID
    :type run_id: int
    :return: Either a JSON-RPC2-compliant success or failure message,
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Cancel Run {params['run_id']}")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.mark_to_cancel()
    except Exception as error:
        return failure(error)
    return success(result)


def submit_run(params, session):
    """Save new run submission in database.  The daemon shall submit to Cromwell after Globus tranfer completes."""
    logger.info(f"User {params['user_id']}: Submit Run {params['run_id']}")
    try:
        run = Run.from_params(session, params)
    except Exception as error:
        return failure(error)
    else:
        return success(run.data.status)


def resubmit_run(params, session):
    """Resubmit a run"""
    logger.info(f"User {params['user_id']}: Re-submit Run {params['run_id']}")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.resubmit()
    except Exception as error:
        return failure(error)
    else:
        return success(result)


def run_task_log(params, session):
    """Retrieve task log from Cromwell metadata"""
    logger.info(f"User {params['user_id']}: Task-log Run {params['run_id']}")
    try:
        task_log = TaskLog(session, params["cromwell_run_id"], logger)
        local_tz = params.get("local_tz", DEFAULT_TZ)
        result = task_log.table(local_tz=local_tz)
    except Exception as error:
        return failure(error)
    else:
        return success(result)


def submit_transfer(params, session):
    """
    Direct the Site to transfer files.  The transfer is added to the queue and will be done later.
    The status of the transfer may be queried using the &transfer_status function below.
    """
    logger.info(f"New transfer {params['transfer_id']}")
    try:
        transfer = Transfer.from_params(session, params)
    except Exception as error:
        logger.debug(f"Error submitting transfer {params['transfer_id']}: {error}")
        return failure(error)
    else:
        result = {"status": transfer.status()}
        return success(result)


def transfer_status(params, session):
    """
    Check the status of a transfer.
    """
    try:
        transfer = Transfer.from_id(session, params["transfer_id"])
    except Exception as error:
        logger.error(f"Transfer {params['transfer_id']} status failed: {error}")
        return failure(error)
    else:
        result = {"status": transfer.status(), "reason": transfer.reason()}
        return success(result)


def cancel_transfer(params, session):
    """
    Check the status of a transfer.
    """
    logger.info(f"Cancel transfer {params['transfer_id']}")
    try:
        transfer = Transfer.from_id(session, params["transfer_id"])
        transfer.cancel()
    except Exception as error:
        logger.debug(f"Error cancelling transfer {params['transfer_id']}: {error}")
        return failure(error)
    else:
        result = {"status": transfer.status()}
        return success(result)


def site_config(params, session):
    """
    Return site parameters to central.
    """
    try:
        result = config.conf.get_site_config()
    except Exception as error:
        logger.error(f"Failed to get site config: {error}")
        return failure(error)
    else:
        return success(result)


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
operations = {
    "server_status": {"function": server_status},
    "queue_wait": {"function": queue_wait},
    "submit_run": {
        "function": submit_run,
        "required_params": [
            "user_id",
            "run_id",
            "submission_id",
            "input_site_id",
        ],
    },
    "resubmit_run": {
        "function": resubmit_run,
        "required_params": ["run_id"],
    },
    "cancel_run": {
        "function": cancel_run,
        "required_params": ["user_id", "run_id"],
    },
    "run_task_log": {
        "function": run_task_log,
        "required_params": ["user_id", "run_id"],
    },
    "submit_transfer": {
        "function": submit_transfer,
        "required_params": [
            "transfer_id",
            "src_base_dir",
            "dest_base_dir",
            "manifest",
        ],
    },
    "transfer_status": {
        "function": transfer_status,
        "required_params": ["transfer_id"],
    },
    "cancel_transfer": {
        "function": cancel_transfer,
        "required_params": ["transfer_id"],
    },
    "site_config": {
        "function": site_config,
        "required_params": [],
    },
}
