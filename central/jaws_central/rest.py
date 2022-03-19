"""
REST endpoints.
"""

import logging
import json
from elasticsearch import Elasticsearch
from flask import abort, request
from jaws_central import config, runs
from jaws_central.runs import Run
from jaws_central.models_fsa import db


logger = logging.getLogger(__package__)


def search_runs(user):
    """Search is used by both user queue and history commands."""
    logger.info(f"User {user}: Search runs")
    result = runs.search_runs(
        db.session,
        user_id=user,
        site_id=request.form.get("site_id", "all").upper(),
        active_only=True if request.form.get("active_only") == "True" else False,
        delta_days=int(request.form.get("delta_days", 0)),
        result=request.form.get("result", "any").lower(),
        all_users=True if request.form.get("all") == "True" else False,
    )
    return result, 200


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
        max_ram_gb = config.conf.get_site(site_id, "max_ram_gb")
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
    result["uploads_dir"] = f'{result["uploads_dir"]}/{user}'
    return result, 200


def submit_run(user):
    """
    Record the run submission in the database, with status as "uploading".

    :param user: current user's ID
    :type user: str
    :return: run_id, upload_task_id
    :rtype: dict
    """
    try:
        current_user = db.session.query(User).get(user)
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, {"error": f"Db error; {e}"})
    try:
        run = Run(
            db.session,
            site_id=request.form.get("site_id", None).upper(),
            submission_id=request.form.get("submission_id"),
            input_site_id=request.form.get("input_site_id", None).upper(),
            input_endpoint=request.form.get("input_endpoint", None),
            output_endpoint=request.form.get("output_endpoint"),
            output_dir=request.form.get("output_dir"),
            wdl_file=request.form.get("wdl_file"),
            json_file=request.form.get("json_file"),
            tag=request.form.get("tag"),
            user_id=user,
            user_email=current_user.email,
        )
    except Exception as error:
        logger.error(f"User {user}: Submit Run failed: {error}")
        abort(500, {"error", f"{error}"})
    else:
        logger.info(f"User {user}: Submit Run {run.id}")
        result = {
            "run_id": run.id,
            "status": run.status,
        }
        return result, 201


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
    logger.info(f"User {user}: Get status of Run {run_id}")
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.status(verbose)
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


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
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.task_status()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


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
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.run_log()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


def task_log(user, run_id):
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
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.task_log()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


def task_summary(user, run_id):
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
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.task_summary()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


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
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.metadata()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


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
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.outputs()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


def get_errors(user, run_id):
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
    run = Run(db.session, run_id=run_id, user_id=user)
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.errors()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 200


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
    run = Run(db.session, run_id=run_id, user_id=user)
    try:
        run = Run(db.session, run_id=run_id, user_id=user)
        result = run.errors()
    except Exception as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    else:
        return result, 201


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
        active_runs = runs.select_runs(db.session, user, active_only=True)
    except SQLAlchemyError as error:
        logger.exception(error)
        abort(500, {"error": f"{error}"})
    results = {}
    for run_id in active_runs:
        result = None
        try:
            run = Run(db.session, run_id=run_id, user_id=user)
            result = run.cancel()
        except Exception as error:
            logger.exception(error)
            result = f"{error}"
        else:
            results[run_id] = result
    return results, 201


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
    logger.info(f"User {user}: Get log of Run {run_id}")
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
