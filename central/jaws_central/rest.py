"""
JAWS-Central REST endpoints.
"""

import logging
from datetime import datetime, timedelta
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config, db, run
from jaws_central.user import User, Run

logger = logging.getLogger(__package__)


def user_queue(user):
    """Return the current user's unfinished runs.

    :param user: current user's ID
    :type user: str
    :return: details about any current runs
    :rtype: list
    """
    logger.info(f"User {user}: Get queue")
    try:
        current_user = User(db.session, user)
        queue = run.queue(db.session, current_user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    return queue, 200


def user_history(user, delta_days=10):
    """Return the current user's recent runs, regardless of status.

    :param user: current user's ID
    :type user: str
    :param delta_days: number of days in which to search
    :type delta_days: int
    :return: details about any recent runs
    :rtype: list
    """
    start_date = datetime.today() - timedelta(int(delta_days))
    logger.info(f"User {user}: Get history, last {delta_days} days")
    try:
        current_user = User(db.session, user)
        history = run.history(db.session, current_user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    return history, 200


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
        abort(404, f'Unknown Site ID; "{site_id}" is not one of our sites')
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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, params)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    logger.debug(f"User {user}: New run {run.id} to {run.model.site_id}")
    result = {
        "run_id": run.id,
        "site_id": run.model.site_id,
        "tag": run.model.tag,
    }
    return result, 201


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    result = {"run_id": run_id, "status": run.status}
    return result, 200


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return run.task_status(), 200


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return run.log(), 200


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return run.task_log(), 200


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return run.metadata(), 200


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return run.errors(), 200


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
        current_user = User(db.session, user)
        run = Run(db.session, current_user, run_id=run_id)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    except RunNotFoundError as error:
        abort(404, f"Run {run_id}: {error}")
    return run.cancel(), 201
