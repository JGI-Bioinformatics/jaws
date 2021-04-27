"""
JAWS-Central REST interface; also see swagger.rest.yml
"""

import logging
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config, db, runs
from jaws_central.runs import Run, RunNotFoundError

logger = logging.getLogger(__package__)


def search_runs(user):
    """Search is used by both user queue and history commands."""
    site_id = request.form.get("site_id", "all").upper()
    active_only = True if request.form.get("active_only") == "True" else False
    delta_days = int(request.form.get("delta_days", 0))
    result = request.form.get("result", "any").lower()
    logger.info(f"User {user}: Search runs")
    search_results = runs.select_runs(user, active_only, delta_days, site_id, result)
    return search_results, 200


def list_sites(user):
    """List all JAWS-Sites.

    :param user: current user's ID
    :type user: str
    :return: available sites and max ram of each
    :rtype: dict
    """
    logger.info(f"User {user}: List sites")
    return config.conf.available_sites(), 200


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
    site_info = config.conf.get_site_info(site_id)
    if site_info is None:
        abort(404, f'Unknown Site ID; "{site_id}" is not one of our sites')
    site_info["uploads_dir"] = f'{site_info["uploads_dir"]}/{user}'
    return site_info, 200


def submit_run(user):
    """
    Record the run submission in the database, with status as "uploading".

    :param user: current user's ID
    :type user: str
    :return: run_id, upload_task_id
    :rtype: dict
    """
    params = {}
    params["site_id"] = request.form.get("site_id", None).upper()
    params["submission_id"] = request.form.get("submission_id")
    params["input_site_id"] = request.form.get("input_site_id", None).upper()
    params["input_endpoint"] = request.form.get("input_endpoint", None)
    params["output_endpoint"] = request.form.get("output_endpoint")
    params["output_dir"] = request.form.get("output_dir")
    params["wdl_file"] = request.form.get("wdl_file")
    params["json_file"] = request.form.get("json_file")
    params["tag"] = request.form.get("tag")

    try:
        run = Run(db.session, user, params)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    logger.debug(f"User {user}: New run {run.id} to {run.model.site_id}")
    run_submission = {
        "run_id": run.id,
        "site_id": run.model.site_id,
        "tag": run.model.tag,
    }
    return run_submission, 201


def run_status(user, run_id):
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
        run = Run(db.session, user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    current_status = {"run_id": run_id, "status": run.status}
    return current_status, 200


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
        run = Run(db.session, user, run_id=run_id)
        task_status = run.task_status()
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return task_status, 200


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
        run = Run(db.session, user, run_id=run_id)
        log = run.log()
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return log, 200


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
        run = Run(db.session, user, run_id=run_id)
        task_log = run.task_log()
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return task_log, 200


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
        run = Run(db.session, user, run_id=run_id)
        cromwell_metadata = run.metadata()
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return cromwell_metadata, 200


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
    try:
        run = Run(db.session, user, run_id=run_id)
        errors_report = run.errors()
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return errors_report, 200


def cancel_run(user, run_id):
    """
    Cancel a run.  It doesn't cancel Globus transfers, just Cromwell runs.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f"User {user}: Cancel Run {run_id}")
    try:
        run = Run(db.session, user, run_id=run_id)
        run.cancel()
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return {"cancel": "OK"}, 201


def cancel_all(user):
    """
    Cancel all of a user's active runs.

    :param user: current user's ID
    :type user: str
    :return: runs that were cancelled
    :rtype: dict
    """
    logger.info(f"User {user}: Cancel-all")
    try:
        runs_cancelled = runs.cancel_all(db.session, user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    return runs_cancelled, 201
