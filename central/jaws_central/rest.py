"""
Analysis (AKA Run) REST endpoints.
"""

import logging
import json
from elasticsearch import Elasticsearch
from flask import abort, request, current_app as app
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config
from jaws_rpc import rpc_index
from jaws_central.runs import Run, select_runs, RunNotFoundError
from jaws_central.users import User
from jaws_central.config import ConfigurationError


logger = logging.getLogger(__package__)


run_active_states = [
    "created",
    "upload queued",
    "uploading",
    "upload inactive",
    "upload complete",
    "submitted",
    "queued",
    "running",
    "succeeded",
    "ready",
    "download queued",
    "downloading",
]


run_pre_cromwell_states = [
    "created",
    "upload queued",
    "uploading",
    "upload inactivte",
    "upload complete",
    "upload stalled",
    "upload failed",
    "ready",
    "submission failed",
]


class RunAccessDeniedError(Exception):
    pass


def rpc_call(user, run_id, method, params={}):
    """This is not a Flask endpoint, but a helper used by several endpoints.
    It checks a user's permission to access a run, perform the specified RPC function,
    and returns result if OK, aborts if error.

    :param user: current user's id
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :param method: the method to execute remotely
    :type method: string
    :param params: parameters for the remote method, depends on method
    :type params: dict
    :return: response in JSON-RPC2 format
    :rtype: dict or list
    """
    try:
        response = _rpc_call(user, run_id, method, params)
    except RunNotFoundError as error:
        abort(404, {"error": f"{error}"})
    except RunAccessDeniedError as error:
        abort(401, {"error": f"{error}"})
    except Exception as error:
        abort(500, {"error": f"{error}"})
    if "error" in response:
        abort(response["error"]["code"], {"error": response["error"]["message"]})
    else:
        return response["result"], 200


def _rpc_call(user, run_id, method, params={}):
    """
    It checks a user's permission to access a run, perform the specified RPC function, and
    returns the response, which may indicate success or failure, to be processed by the caller.

    :param user: current user's id
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :param method: the method to execute remotely
    :type method: string
    :param params: parameters for the remote method, depends on method
    :type params: dict
    :return: response in JSON-RPC2 format
    :rtype: dict or list
    """
    try:
        run = Run.from_id(app.session, run_id)
    except RunNotFoundError:
        raise
    except Exception as error:
        logger.error(f"Unable to init Run {run_id}: {error}")
        raise
    if run.data.id != user and not _is_admin(user):
        raise RunAccessDeniedError(
            "Access denied; you cannot access to another user's workflow"
        )
    a_site_rpc_client = rpc_index.rpc_index.get_client(run.data.compute_site_id)
    params["user_id"] = user
    params["run_id"] = run_id
    params["cromwell_run_id"] = run.data.cromwell_run_id
    logger.info(f"User {user} RPC {method} params {params}")
    try:
        response = a_site_rpc_client.request(method, params)
    except Exception as error:
        logger.error(f"RPC {method} failed: {error}")
        raise
    return response


def _get_user(user):
    """
    Retrieve current user info.
    :param user: Current user's ID
    :type user: str
    :return: User object
    :rtype: users.User
    """
    try:
        current_user = User.from_id(app.session, user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Error retrieving user record: {error}"})
    return current_user


def _is_admin(user):
    """
    Check if current user is an adinistrator.
    :param user: Current user's ID
    :type user: str
    :return: True if user is an admin
    :rtype: bool
    """
    current_user = _get_user(user)
    return True if current_user.data.is_admin else False


def _user_group(user):
    """
    Check if current user belongs to a user-group.
    :param user: Current user's ID
    :type user: str
    :return: user-group (may be None)
    :rtype: str
    """
    current_user = _get_user(user)
    return current_user.data.user_group


def search_runs(user, verbose=False):
    """Search is used by both user queue and history commands."""
    site_id = request.form.get("site_id", "all").upper()
    active_only = True if request.form.get("active_only") == "True" else False
    delta_days = int(request.form.get("delta_days", 0))
    result = request.form.get("result", "any").lower()
    all_users = True if request.form.get("all") == "True" else False
    logger.info(f"User {user}: Search runs")
    matches = select_runs(
        app.session,
        user_id=user,
        active_only=active_only,
        delta_days=delta_days,
        site_id=site_id,
        result=result,
        all_users=all_users,
    )
    runs = []
    for run in matches:
        runs.append(run.info(verbose))
    return runs, 200


def list_sites(user):
    """List all JAWS-Sites.

    :param user: current user's ID
    :type user: str
    :return: list of valid Site IDs
    :rtype: list
    """
    logger.info(f"User {user}: List sites")
    result = []
    for site_id in config.conf.sites.keys():
        max_ram_gb = config.conf.get_site_param(site_id, "max_ram_gb")
        record = {
            "site_id": site_id,
            "max_ram_gb": max_ram_gb,
        }
        result.append(record)
    return result, 200


def get_site(user, site_id):
    """Get parameters of a Site, required to submit a run.

    :param user: current user's ID
    :type user: str
    :param site_id: a JAWS-Site's ID
    :type site_id: str
    :return: globus endpoint id and uploads path
    :rtype: dict
    """
    logger.debug(f"User {user}: Get info for site {site_id}")
    result = config.conf.get_site_info(site_id)
    if result is None:
        abort(404, {"error": f'Unknown Site ID; "{site_id}" is not one of our sites'})
    result["inputs_dir"] = f'{result["inputs_dir"]}/{user}'
    return result, 200


def submit_run(user):
    """
    Insert the run submission in the database and return the run_id (primary key).

    :param user: current user's ID
    :type user: str
    :return: run metadata including run_id
    :rtype: dict
    """
    input_site_id = request.form.get("input_site_id", None).upper()
    compute_site_id = request.form.get("compute_site_id", None).upper()
    submission_id = request.form.get("submission_id")
    caching = False if request.form.get("caching") == "False" else True
    max_ram_gb = int(request.form.get("max_ram_gb"))
    wdl_file = request.form.get("wdl_file")
    json_file = request.form.get("json_file")
    tag = request.form.get("tag")
    webhook = request.form.get("webhook", None)
    manifest_json = request.form.get("manifest")
    logger.info(
        f"User {user}: New run submission {submission_id} from {input_site_id} to {compute_site_id}"
    )

    # check if valid compute site was requested
    try:
        compute_site_config = config.conf.get_site(compute_site_id)
    except ConfigurationError as error:
        logger.warning(
            f"User requested unknown compute-site, {compute_site_id}: {error}"
        )
        abort(
            404,
            {"error": f'Unknown Site ID; "{compute_site_id}" is not one of our sites'},
        )
    if "user_group" in compute_site_config and compute_site_config["user_group"]:
        # a jaws-site may optionally be restricted to members of a user-group
        current_user = _get_user(app.session, user)
        if (
            current_user.data.user_group == compute_site_config["user_group"]
            or current_user.data.is_admin
        ):
            pass
        else:
            abort(
                401,
                {
                    "error": f"Access denied; use of jaws-{compute_site_id} is restricted"
                },
            )

    # check if requested compute site can process this WDL
    if (
        "max_ram_gb" in compute_site_config
        and int(compute_site_config["max_ram_gb"]) < max_ram_gb
    ):
        msg = f'The requested site {compute_site_id} has only {compute_site_config["max_ram_gb"]} GB which is less than the {max_ram_gb} GB required by this workflow.'  # noqa
        abort(406, {"error": msg})

    # insert into db and the daemon will pick it up and initiate the file transfer
    params = {
        "user_id": user,
        "submission_id": submission_id,
        "status": "created",
        "max_ram_gb": max_ram_gb,
        "caching": caching,
        "input_site_id": input_site_id,
        "compute_site_id": compute_site_id,
        "wdl_file": wdl_file,
        "json_file": json_file,
        "tag": tag,
        "manifest_json": manifest_json,
        "webhook": webhook,
    }
    try:
        run = Run.from_params(app.session, params)
    except Exception as error:
        logger.exception(f"Error inserting Run: {error}")
        abort(500, {"error": f"Error inserting Run into db: {error}"})
    logger.debug(f"User {user}: New run {run.data.id}")

    # return run_id
    result = {
        "run_id": run.data.id,
    }
    logger.info(f"User {user}: New run {run.data.id}")
    return result, 201


def _get_run(user, run_id):
    """Return a Run object if found and accessible by current user, else abort.

    :param run_id: Run primary key
    :type run_id: str
    :return: the run ORM object for the record
    :rtype: sqlalchemy.model
    """
    try:
        run = Run.from_id(app.session, run_id)
    except RunNotFoundError as error:
        logger.error(f"Run {run_id} not found: {error}")
        abort(404, {"error": "Run not found; please check your run_id"})
    except Exception as error:
        logger.error(error)
        abort(500, {"error": f"Error retrieving Run {run_id}; {error}"})
    if run.data.user_id != user and not _is_admin(user):
        abort(401, {"error": "Access denied; you are not the owner of that Run."})
    return run


def _abort_if_pre_cromwell(run):
    """Returns if run was submitted to Cromwell, aborts otherwise.

    :param run: SQLAlchemy Model Run object
    :type param: sqlalchemy.model
    :return: True if pre-cromwell submission; false otherwise.
    :rtype: boolean
    """
    if run.status in run_pre_cromwell_states:
        abort(
            404,
            {
                "error": "No data available as the Run hasn't been submitted to Cromwell yet."
            },
        )


def run_status(user, run_id, verbose=False):
    """
    Retrieve the current status of a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: The status of the run, if found; abort otherwise
    :rtype: dict
    """
    run = _get_run(user, run_id)
    logger.info(f"User {user}: Get status of Run {run.data.id}")
    return run.info(verbose), 200


def run_status_complete(user, run_id):
    """
    Retrieve the current status of a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: The status of the run, if found; abort otherwise
    :rtype: dict
    """
    return run_status(user, run_id, True)


def task_status(user, run_id):
    """
    Retrieve the current status of each Task in a Run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: The status of each task in a run.
    :rtype: dict
    """
    logger.info(f"User {user}: Get task-status of Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return rpc_call(user, run_id, "get_task_status")


def run_log(user: str, run_id: int):
    """
    Retrieve complete log of a Run's state transitions.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a Run
    :type run_id: int
    :return: Table of log entries
    :rtype: list
    """
    logger.info(f"User {user}: Get log of Run {run_id}")
    run = _get_run(user, run_id)
    try:
        log = run.log()
    except Exception as error:
        logger.exception(f"Error retrieving log of run {run_id}: {error}")
        abort(500, {"error": str(error)})
    else:
        return log, 200


def run_task_log(user, run_id):
    """
    Retrieve log of all task state transitions.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: The complete log of all task state transitions.
    :rtype: dict
    """
    logger.info(f"User {user}: Get task-log for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return rpc_call(user, run_id, "run_task_log")


def run_task_summary(user, run_id):
    """
    Retrieve summary of all task state transitions.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: The complete summary of all task state transitions.
    :rtype: dict
    """
    logger.info(f"User {user}: Get task-log for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return rpc_call(user, run_id, "run_task_summary")


def run_metadata(user, run_id):
    """
    Retrieve the metadata of a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell metadata for the run, if found; abort otherwise
    :rtype: dict
    """
    logger.info(f"User {user}: Get metadata for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return rpc_call(user, run_id, "run_metadata")


def run_outputs(user, run_id):
    """
    Retrieve the outputs of a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell outputs for the run, if found; abort otherwise
    :rtype: dict
    """
    logger.info(f"User {user}: Get outputs for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    result = rpc_call(user, run_id, "run_outputs")
    return result


def run_outfiles(user, run_id):
    """
    Retrieve the output files of a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: List of output files
    :rtype: list
    """
    logger.info(f"User {user}: Get outfiles for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    result = rpc_call(user, run_id, "run_outfiles")
    return result


def run_workflow_root(user, run_id):
    """
    Retrieve the workflow_root for the Run.  It refers to the path at the compute-site.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Unmodified root folder of the Run
    :rtype: str
    """
    logger.info(f"User {user}: Get outfiles for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    result = rpc_call(user, run_id, "run_workflow_root")
    return result


def run_errors(user, run_id):
    """
    Retrieve error messages and stderr for failed tasks.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell error messages and stderr for all failed tasks.
    :rtype: str
    """
    logger.info(f"User {user}: Get errors for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return rpc_call(user, run_id, "run_errors")


def run_running_tasks(user, run_id):
    """
    Retrieve running-tasks report.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell running tasks report.
    :rtype: str
    """
    logger.info(f"User {user}: Get running-tasks report for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return rpc_call(user, run_id, "run_running_tasks")


def cancel_run(user, run_id):
    """
    Cancel a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f"User {user}: Cancel Run {run_id}")

    run = _get_run(user, run_id)
    try:
        run.cancel()
    except Exception as error:
        abort(400, {"error": str(error)})
    else:
        return {run_id: True}, 201


def cancel_all(user):
    """
    Cancel all of a user's active runs.

    :param user: current user's ID
    :type user: str
    :return: run ids and results
    :rtype: dict
    """
    logger.info(f"User {user}: Cancel-all")
    try:
        active_runs = select_runs(app.session, user_id=user, active_only=True)
    except Exception as error:
        logger.error(error)
        abort(500, {"error": f"Search runs error: {error}"})
    cancelled = {}
    for run in active_runs:
        try:
            run.cancel()
        except Exception as error:
            logger.warn(f"Error cancelling run {run.data.id}: {error}")
            cancelled[run.data.id] = False
        else:
            cancelled[run.data.id] = True
    return cancelled, 201


def _search_elastic_search(host, port, api_key, index, query, aggregations=None):
    """
    Search Elastic Search (ES) DB

    :param elastic_client: ES Client
    :type elastic_client: obj
    :param query: ES top-level query
    :type query: dict
    :param aggregations: ES bucket aggregations
    :type aggregations: dict
    :return: ES search response
    :rtype: dict
    """
    try:
        elastic_client = Elasticsearch([f"http://{host}:{port}"], api_key=api_key)
        response = elastic_client.search(
            index=index, query=query, aggregations=aggregations, size=10000
        )
    except Elasticsearch.AuthorizationException as error:
        logger.error(error)
        abort(403, {"error": f"Not authorized; {error}"})
    except Elasticsearch.AuthenticationException as error:
        logger.error(error)
        abort(401, {"error": f"Invalid/Missing credentials; {error}"})
    except Elasticsearch.NotFoundError as error:
        logger.error(error)
        abort(404, {"error": f"Not found; {error}"})
    except Exception as error:
        logger.error(error)
        abort(500, {"error": f"Elastic-Search error; {error}"})
    return response


def get_performance_metrics(user, run_id):
    """
    Query ES to get performance metrics for user's run

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: performance metrics
    :rtype: dict
    """
    run = _get_run(user, run_id)
    logger.info(f"User {user}: Get log of Run {run.data.id}")
    db_conf = config.conf.get_section("ELASTIC_SEARCH")
    pm_conf = config.conf.get_section("PERFORMANCE_METRICS")
    response = _search_elastic_search(
        host=db_conf["host"],
        port=db_conf["port"],
        api_key=db_conf["api_key"],
        index=pm_conf["index"],
        query={"match": {"jaws_run_id": int(run_id)}},
    )
    metrics = []
    try:
        for hit in response["hits"]["hits"]:
            metrics.append(hit["_source"])
    except Exception as error:
        status = 404
        if "error" in response and "status" in response:
            status = response["status"]
            error += json.dumps(response["error"], indent=2)
        logger.error(error)
        abort(status, {"error": f"{error}"})
    return metrics
