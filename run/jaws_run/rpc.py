"""
JAWS Run Service RPC methods.
"""

import logging
from jaws_run.config import conf
from jaws_run.api.run import (
    Run,
    RunNotFoundError,
    RunAlreadyExistsError,
    DatabaseError,
    CromwellError,
    TaskServiceError,
)
from jaws_run.api.cromwell import Cromwell
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)
cromwell = Cromwell(conf.get("CROMWELL", "url"))


def get_cromwell_status(params):
    """
    Return the current status of the Cromwell service.

    :param params: None
    :type params: dict
    :return: JSON-RPC2 successful response with no content
    :rtype: dict
    """
    logger.debug("Check Cromwell status")
    try:
        _ = cromwell.status()
    except Exception:
        return failure(500, "Cromwell offline")
    return success()


def submit(params):
    """
    Submit new run for execution.

    :param params: run_id, user_id, submission_id, output_endpoint, output_dir
    :type params: dict
    :return: JSON-RPC2 response to no content
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        _ = Run(run_id, params)
    except RunAlreadyExistsError:
        return failure(404, f"Run {run_id} already exists")
    except DatabaseError as error:
        return failure(500, f"Run db offline: {error}")
    except ValueError as error:
        return failure(401, f"Invalid Run parameters: {error}")
    return success()


# TODO THIS WAS get_status, NEED TO UPDATE jaws-central
def get_runs(params):
    """
    Get the run info, including current status, for a list of Runs.
    Runs may be specified by run_id and/or cromwell_run_id.
    This is used by jaws-central to get info for a single run or status for a list of runs.

    :param params: List of run_ids
    :type params: dict
    :return: JSON-RPC2 response with result of run_id:status
    :rtype: dict
    """
    results = {}
    for run_id in params["run_ids"]:
        try:
            run = Run(run_id=uun_id)
            result = run.get_status()
        except RunNotFoundError:
            return failure(404, f"Run {run_id} not found")
        except Exception as error:
            return failure(500, f"{error}")
        results.append(result)
    return success(results)


def get_log(params):
    """
    Get a log of a state transitions of a Run.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with table of Run state transitions.
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
        result = run.get_log()
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)


def get_metadata(params):
    """
    Retrieve the Cromwell metadata of a run and all it's subworkflows.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with dict of cromwell_run_id:Cromwell metadata JSON.
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
        result = run.get_metadata()
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)


def cancel(params):
    """
    Cancel a run.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with no content
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
        run.cancel()
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except Exception as error:
        return failure(500, f"{error}")
    return success()


def get_errors(params):
    """
    Extract error messages from Cromwell metadata and include stderr for failed Tasks.

    :param params: run_id
    :type params: dict
    :return: JSON-RPC2 response with dict of failed tasks: err msg and stderr
    :rtype: dict
    """
    run_id = params["run_id"]
    try:
        run = Run(run_id)
        result = run.get_errors()
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)


def get_tasks(params):
    """
    Get all task_ids associated with a Run.

    :param params: must include "cromwell_run_id"
    :type params: dict
    :return: "result" shall be list of run_ids 
    :rtype: dict
    """
    try:
        run = Run(run_id)
        result = run.get_tasks()  # TODO add to API
    except RunNotFoundError:
        return failure(404, f"Run {run_id} not found")
    except Exception as error:
        return failure(500, f"{error}")
    return success(result)

# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
rpc_methods = {
    "get_cromwell_status": {"method": get_cromwell_status},
    "submit": {
        "method": submit,
        "required_params": [
            "user_id",
            "run_id",
            "email",
            "transfer_refresh_token",
            "submission_id",
            "upload_task_id",
            "output_endpoint",
            "output_dir",
        ],
    },
    "get_runs": {"method": get_runs, "required_params": ["run_ids"]},
    "get_log": {"method": get_log, "required_params": ["run_id"]},
    "get_metadata": {"method": get_metadata, "required_params": ["run_id"]},
    "cancel": {"method": cancel, "required_params": ["run_id"]},
    "get_errors": {"method": get_errors, "required_params": ["run_id"]},
    "get_tasks": {"method": get_tasks, "required_params": ["cromwell_run_id"]},
}
