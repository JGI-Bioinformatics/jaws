"""
Analysis (AKA Run) REST endpoints.
"""

import logging
import json
from datetime import datetime, timedelta
from elasticsearch import Elasticsearch
from flask import abort, request, current_app as app
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config, jaws_constants
from jaws_rpc import rpc_index
from jaws_central.models import Run, User, Run_Log
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


class RunNotFoundError(Exception):
    pass


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
        run = app.session.query(Run).get(run_id)
    except SQLAlchemyError as e:
        logger.error(e)
        raise
    if not run:
        raise RunNotFoundError("Run not found; please check your run_id")
    if run.user_id != user and not _is_admin(user):
        raise RunAccessDeniedError(
            "Access denied; you cannot access to another user's workflow"
        )
    a_site_rpc_client = rpc_index.rpc_index.get_client(run.compute_site_id)
    params["user_id"] = user
    params["run_id"] = run_id
    params["cromwell_run_id"] = run.cromwell_run_id
    logger.info(f"User {user} RPC {method} params {params}")
    try:
        response = a_site_rpc_client.request(method, params)
    except Exception as error:
        logger.error(f"RPC {method} failed: {error}")
        raise
    return response


def _is_admin(user):
    """
    Check if current user is an adinistrator.
    :param user: Current user's ID
    :type user: str
    :return: True if user is an admin
    :rtype: bool
    """
    try:
        current_user = app.session.query(User).get(user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Db error; {error}"})
    return True if current_user.is_admin else False


def _user_group(user):
    """
    Check if current user belongs to a user-group.
    :param user: Current user's ID
    :type user: str
    :return: user-group (may be None)
    :rtype: str
    """
    try:
        current_user = app.session.query(User).get(user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Db error; {error}"})
    return current_user.user_group


def search_runs(user):
    """Search is used by both user queue and history commands."""
    site_id = request.form.get("site_id", "all").upper()
    active_only = True if request.form.get("active_only") == "True" else False
    delta_days = int(request.form.get("delta_days", 0))
    result = request.form.get("result", "any").lower()
    all_users = True if request.form.get("all") == "True" else False
    logger.info(f"User {user}: Search runs")
    is_admin = _is_admin(user)
    rows = _select_runs(
        user,
        active_only=active_only,
        delta_days=delta_days,
        site_id=site_id,
        result=result,
        all_users=all_users,
    )
    runs = []
    for row in rows:
        runs.append(_run_info(row, is_admin))
    return runs, 200


def _select_runs(user: str, **kwargs):
    """Select runs from db.

    :param user: current user's ID
    :type user: str
    :return: Runs matching search criteria
    :rtype: list
    """
    query = app.session.query(Run)
    if "all_users" in kwargs and kwargs["all_users"] is True:
        pass
    else:
        query = query.filter(Run.user_id == user)
    if "active_only" in kwargs and kwargs["active_only"] is True:
        query = query.filter(Run.status.in_(run_active_states))
    if "site_id" in kwargs:
        site_id = kwargs["site_id"].upper()
        if site_id != "ALL":
            query = query.filter(Run.compute_site_id == site_id)
    if "delta_days" in kwargs:
        delta_days = int(kwargs["delta_days"])
        if delta_days > 0:
            start_date = datetime.today() - timedelta(delta_days)
            query = query.filter(Run.submitted >= start_date)
    if "result" in kwargs:
        result = kwargs["result"].lower()
        if result != "any":
            query = query.filter(Run.result == result)
    return query.all()


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
    if (
        "user_group" in compute_site_config
        and compute_site_config["user_group"]
    ):
        # a jaws-site may optionally be restricted to members of a user-group
        user_group = _user_group(user)
        if user_group == compute_site_config["user_group"] or _is_admin(user):
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
    # We aren't using the Run class here because it isn't compatible with the
    # flask-sqlalchemy base class
    # TODO: stop using flask-sqlalchemy and use runs.Run.from_params() instead
    run = Run(
        user_id=user,
        submission_id=submission_id,
        status="created",
        max_ram_gb=max_ram_gb,
        caching=caching,
        input_site_id=input_site_id,
        compute_site_id=compute_site_id,
        wdl_file=wdl_file,
        json_file=json_file,
        tag=tag,
        manifest_json=manifest_json,
        webhook=webhook,
    )
    try:
        app.session.add(run)
    except Exception as error:
        app.session.rollback()
        logger.exception(f"Error inserting Run: {error}")
        abort(500, {"error": f"Error inserting Run into db: {error}"})
    try:
        app.session.commit()
    except Exception as error:
        app.session.rollback()
        err_msg = f"Unable to insert new run in db: {error}"
        logger.exception(err_msg)
        abort(500, {"error": err_msg})
    logger.debug(f"User {user}: New run {run.id}")

    # return run_id
    result = {
        "run_id": run.id,
    }
    logger.info(f"User {user}: New run {run.id}")
    return result, 201


def _update_run_status(run, new_status, reason=None):
    """Update run table and insert run_logs entry."""
    status_from = run.status
    run.status = new_status
    if new_status == "cancelled":
        run.result = "cancelled"
    try:
        app.session.commit()
    except Exception as error:
        app.session.rollback()
        logger.exception(f"Error updating run status in db: {error}")
    log = Run_Log(
        run_id=run.id,
        status_from=status_from,
        status_to=run.status,
        timestamp=run.updated,
        reason=reason,
    )
    try:
        app.session.add(log)
        app.session.commit()
    except Exception as error:
        app.session.rollback()
        logger.exception(f"Error insert run log entry: {error}")


def _get_run(user, run_id):
    """Return a Run object if found and accessible by current user, else abort.

    :param run_id: Run primary key
    :type run_id: str
    :return: the run ORM object for the record
    :rtype: sqlalchemy.model
    """
    try:
        run = app.session.query(Run).get(run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Db error; {error}"})
    if not run:
        abort(404, {"error": "Run not found; please check your run_id"})
    if run.user_id != user and not _is_admin(user):
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
    logger.info(f"User {user}: Get status of Run {run.id}")
    is_admin = _is_admin(user)
    info = _run_info(run, is_admin, verbose)
    return info, 200


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
    run = _get_run(user, run_id)
    logger.info(f"User {user}: Get log of Run {run.id}")
    try:
        query = (
            app.session.query(Run_Log)
            .filter_by(run_id=run_id)
            .order_by(Run_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting from run_logs: {error}")
        abort(500, {"error": f"Db error; {error}"})
    table = []
    for log in query:
        reason = log.reason if log.reason else ""
        row = [
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            reason,
        ]
        table.append(row)
    return table, 200


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
    status = run.status

    if status == "cancelled":
        abort(400, {"error": "That Run had already been cancelled"})
    elif status in ["download complete", "email sent", "done"]:
        abort(400, {"error": "It is too late to cancel."})
    cancelled = _cancel_run(user, run)
    if cancelled is True:
        return {run_id: cancelled}, 201
    else:
        abort(400, {"error": f"Run {run_id} could not be cancelled"})


def _cancel_run(user, run, reason="Cancelled by user"):
    """
    Cancel a Run.

    :param run: Run SqlAlchemy ORM object
    :type run: obj
    """
    orig_status = run.status
    if run.status in ["cancelled", "download complete", "email sent", "done"]:
        # too late to cancel
        return False
    else:
        # central must be updated so the central-run-daemon will not process it
        try:
            _update_run_status(run, "cancelled", reason)
        except Exception as error:
            logger.error(f"Cancel failed to update Run {run.id} status: {error}")
            return False
    if orig_status in ["submitted", "queued", "running"]:
        # the compute jaws-site should also receive the cancel instruction
        try:
            _rpc_call(user, run.id, "cancel_run")
        except Exception as error:
            logger.error(f"RPC cancel Run {run.id} failed: {error}")
            return False
        else:
            logger.debug(f"RPC cancel Run {run.id} successful")
            return True
    else:
        return True


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
        active_runs = _select_runs(user, active_only=True)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Db error; {error}"})
    cancelled = {}
    for run in active_runs:
        cancelled[run.id] = _cancel_run(user, run)
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
    logger.info(f"User {user}: Get log of Run {run.id}")
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


def _run_info(run, is_admin=False, verbose=False):
    """
    Return dictionary of run info.
    Run run cannot be changed by altering the returned dict.
    :param verbose: True if more fields desired else fewer.
    :type verbose: bool
    :return: selected fields
    :rtype: dict
    """
    info = {
        "id": run.id,
        "result": run.result,
        "status": run.status,
        "status_detail": jaws_constants.run_status_msg.get(run.status, ""),
        "compute_site_id": run.compute_site_id,
        "submitted": run.submitted.strftime("%Y-%m-%d %H:%M:%S"),
        "updated": run.updated.strftime("%Y-%m-%d %H:%M:%S"),
        "tag": run.tag,
        "wdl_file": run.wdl_file,
        "json_file": run.json_file,
    }
    if verbose or is_admin:
        more_info = {
            "cromwell_run_id": run.cromwell_run_id,
            "input_site_id": run.input_site_id,
            "upload_id": run.upload_id,
            "submission_id": run.submission_id,
            "download_id": run.download_id,
            "user_id": run.user_id,
        }
        info.update(more_info)
    return info
