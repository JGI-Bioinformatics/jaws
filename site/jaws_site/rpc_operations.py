import logging
from jaws_rpc.responses import success, failure
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.runs import Run, RunNotFoundError
from jaws_site.transfers import Transfer


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


def run_metadata(params, session):
    """Retrieve the metadata of a run.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: The Cromwell metadata for the specified run.
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Metadata Run {params['run_id']}")
    try:
        run = Run.from_id(session, params["run_id"])
        metadata = run.metadata()
        result = metadata.data
    except Exception as error:
        return failure(error)
    return success(result)


def run_outputs(params, session):
    """Retrieve the outputs-json of a run.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: The Cromwell outputs for the specified run.
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Outputs Run {params['run_id']}")
    relpath = False if "relpath" in params and params["relpath"] is False else True
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.outputs(relpath)
    except Exception as error:
        return failure(error)
    return success(result)


def run_outfiles(params, session):
    """Retrieve list of output files of a run.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: The output files for a run
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Outfiles Run {params['run_id']}")
    relpath = False if "relpath" in params and params["relpath"] is False else True
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.outfiles(complete=False, relpath=relpath)
    except Exception as error:
        return failure(error)
    return success(result)


def run_workflow_root(params, session):
    """Retrieve the root dir of the workflow

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: The output files for a run
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: workflowRoot Run {params['run_id']}")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.workflow_root()
    except Exception as error:
        return failure(error)
    return success(result)


def run_manifest(params, session):
    """Retrieve list of output files of a Run.

    :param run_id: JAWS Run ID
    :type params: dict
    :return: The workflow_root and list of output files
    :rtype: list
    """
    logger.info(f"Output manifest for Run {params['run_id']}")
    complete = True if "complete" in params and params["complete"] is True else False
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.output_manifest(complete)
    except Exception as error:
        return failure(error)
    return success(result)


def cancel_run(params, session):
    """Mark a Run to be cancelled.  It will be cancelled by the Run daemon later (asynchronously).

    :param cromwell_run_id: The Cromwell run ID
    :type cromwell_run_id: str
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


def run_errors(params, session):
    """Retrieve error report which is generated from Cromwell metadata.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: errors report
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Run {params['run_id']} errors report")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.errors()
    except RunNotFoundError:
        return success(
            "The Run hasn't yet been submitted to Cromwell so there is no errors report"
        )
    except Exception as error:
        return failure(error)
    return success(result)


def run_running_tasks(params, session):
    """Retrieve running-tasks report which is generated from Cromwell metadata.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: running tasks report
    :rtype: dict
    """
    logger.info(f"User {params['user_id']}: Run {params['run_id']} running-tasks report")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.running_tasks()
    except RunNotFoundError:
        return success(
            "The Run hasn't yet been submitted to Cromwell so there is no running-tasks report"
        )
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


def run_task_log(params, session):
    """Retrieve task log from Cromwell metadata"""
    logger.info(f"User {params['user_id']}: Task-log Run {params['run_id']}")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.task_log()
    except Exception as error:
        return failure(error)
    else:
        return success(result)


def run_task_summary(params, session):
    """Retrieve task summary from Cromwell metadata"""
    logger.info(f"User {params['user_id']}: Task-summary Run {params['run_id']}")
    try:
        run = Run.from_id(session, params["run_id"])
        result = run.task_summary()
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
        result = {
            "status": transfer.status(),
            "reason": transfer.reason()
        }
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
    "run_metadata": {
        "function": run_metadata,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_outputs": {
        "function": run_outputs,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_outfiles": {
        "function": run_outfiles,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_workflow_root": {
        "function": run_workflow_root,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_manifest": {
        "function": run_manifest,
        "required_params": ["run_id"],
    },
    "cancel_run": {
        "function": cancel_run,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_task_log": {
        "function": run_task_log,
        "required_params": ["user_id", "run_id"],
    },
    "run_task_summary": {
        "function": run_task_summary,
        "required_params": ["user_id", "run_id"],
    },
    "run_errors": {
        "function": run_errors,
        "required_params": ["user_id", "cromwell_run_id"],
    },
    "run_running_tasks": {
        "function": run_running_tasks,
        "required_params": ["user_id", "cromwell_run_id"],
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
}
