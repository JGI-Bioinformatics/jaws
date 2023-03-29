"""
This module provides a collection of functions called by rpc_operations.py.
"""

import os
import re
import subprocess

def check_queue_wait(logger):
    """
    Check the queue wait times for all the sites by using the "sbatch --test-only" command which returns an estimatated
    start time for a job. The commands used here are asking for the same resources as condor asks for when creating slurm pools.
    Each site has a combination of small, medium, large, xlarge pools.

    :param logger: The logger passed from rpc_operations.py
    :type params: object
    :return: estimated start times for sm, med, lg & xlg pools in the format: 2023-03-29T10:40:28.
    :rtype: json
    """

    logger.debug("Getting path of queue_wait.sh shim")
    os.getenv('JAWS_BIN_DIR')
    queue_wait_shim = "$JAWS_BIN_DIR/queue-wait.sh"

    try:
        results = subprocess.run([queue_wait_shim], capture_output=True, check=True)
    except CalledProcessError as error:
        logger.error(error)
        return failure(error)

    x = re.search(r"to\sstart\sat\s(\S+)", results.stderr)
    print(x.group(1))

    return True



def server_status(params, session):
    """Return the current status of the Cromwell server.

    :return: Either a success- or failure-formatted JSON-RPC2 response,
    if Cromwell up or not.
    :rtype: dict
    """
    logger.info("Check server status")
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    status = {}
    try:
        status["cromwell"] = cromwell.status()
        status["slurm"] = slurm.status(logger)
        logger.info(status)
    except Exception as error:
        return failure(error)
    return success(status)


def output_manifest(params, session):
    """Retrieve a Run's output manifest (files to return to user).

    :param run_id: JAWS Run ID
    :type params: dict
    :return: The workflow_root and list of output files
    :rtype: list
    """
    logger.info(f"Outfiles for Run {params['run_id']}")
    complete = True if "complete" in params and params["complete"] is True else False
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.output_manifest(complete=complete)
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
        result = task_log.table()
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
    "output_manifest": {
        "function": output_manifest,
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
