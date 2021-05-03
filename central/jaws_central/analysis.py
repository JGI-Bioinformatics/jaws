"""
Analysis (AKA Run) REST endpoints.
"""

import logging
from datetime import datetime, timedelta
from flask import abort, request
from sqlalchemy import and_
from sqlalchemy.exc import SQLAlchemyError
import globus_sdk

import jaws_central.globus

from jaws_central import config
from jaws_central import jaws_constants
from jaws_rpc import rpc_index
from jaws_central.models_fsa import db, Run, User, Run_Log


logger = logging.getLogger(__package__)


run_active_states = [
    "created",
    "uploading",
    "upload inactive",
    "upload complete",
    "submitted",
    "queued",
    "running",
    "succeeded",
    "ready",
    "downloading",
]


run_pre_cromwell_states = [
    "created",
    "uploading",
    "upload inactivte",
    "upload complete",
    "upload stalled",
    "upload failed",
]


def _rpc_call(user, run_id, method, params={}):
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
        run = db.session.query(Run).get(run_id)
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, {"error": f"Db error; {e}"})
    if not run:
        abort(404, {"error": "Run not found; please check your run_id"})
    if run.user_id != user and not _is_admin(user):
        abort(401, {"error": "Access denied; you cannot access to another user's workflow"})
    a_site_rpc_client = rpc_index.rpc_index.get_client(run.site_id)
    params["user_id"] = user
    params["run_id"] = run_id
    params["cromwell_run_id"] = run.cromwell_run_id
    logger.info(f"User {user} RPC {method} params {params}")
    try:
        result = a_site_rpc_client.request(method, params)
    except Exception as error:
        logger.exception(f"RPC {method} failed: {error}")
    if "error" in result:
        abort(result["error"]["code"], {"error": result["error"]["message"]})
    return result["result"], 200


def _is_admin(user):
    """
    Check if current user is an adinistrator.
    :param user: Current user's ID
    :type user: str
    :return: True if user is an admin
    :rtype: bool
    """
    try:
        current_user = db.session.query(User).get(user)
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Db error; {error}"})
    return True if current_user.is_admin else False


def _run_info(run, complete=False):
    """
    Given a SQLAlchemy model for a Run, create a dict with the desired fields.
    :param run: Run object
    :type run: model
    :param complete: True if all fields desired
    :type complete: bool
    :return: selected fields
    :rtype: dict
    """
    info = {}
    if complete:
        info = {
            "id": run.id,
            "submission_id": run.submission_id,
            "cromwell_run_id": run.cromwell_run_id,
            "result": run.result,
            "status": run.status,
            "status_detail": jaws_constants.run_status_msg.get(run.status, ""),
            "site_id": run.site_id,
            "submitted": run.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": run.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "input_site_id": run.input_site_id,
            "input_endpoint": run.input_endpoint,
            "upload_task_id": run.upload_task_id,
            "output_endpoint": run.output_endpoint,
            "output_dir": run.output_dir,
            "download_task_id": run.download_task_id,
            "user_id": run.user_id,
            "tag": run.tag,
            "wdl_file": run.wdl_file,
            "json_file": run.json_file,
        }
    else:
        info = {
            "id": run.id,
            "result": run.result,
            "status": run.status,
            "status_detail": jaws_constants.run_status_msg.get(run.status, ""),
            "site_id": run.site_id,
            "submitted": run.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": run.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "input_site_id": run.input_site_id,
            "tag": run.tag,
            "wdl_file": run.wdl_file,
            "json_file": run.json_file,
        }
    return info


def search_runs(user):
    """Search is used by both user queue and history commands."""
    site_id = request.form.get("site_id", "all").upper()
    active_only = True if request.form.get("active_only") == "True" else False
    delta_days = int(request.form.get("delta_days", 0))
    result = request.form.get("result", "any").lower()
    logger.info(f"User {user}: Search runs")
    is_admin = _is_admin(user)
    """Now passes is_admin to _select_runs"""
    rows = _select_runs(user, active_only, delta_days, site_id, result, is_admin)
    runs = []
    for run in rows:
        runs.append(_run_info(run, is_admin))
    return runs, 200


def _select_runs(user: str, active_only: bool, delta_days: int, site_id: str, result: str, is_admin: bool):
    """Select runs from db.

    :param user: current user's ID
    :type user: str
    :param active_only: Select only active runs
    :type active_only: bool
    :param delta_days: Limit to this many recent days
    :type delta_days: int
    :param site_id: Limit to this compute-site
    :type site_id: str
    :param result: Limit to this result
    :type result: str
    :return: Runs matching search criteria
    :rtype: list
    """
    if not is_admin:
        query = db.session.query(Run).filter(Run.user_id == user)
    else:
        query = db.session.query(Run)
    if active_only:
        query = query.filter(Run.status.in_(run_active_states))
    if site_id != "ALL":
        query = query.filter(Run.site_id == site_id)
    if delta_days > 0:
        start_date = datetime.today() - timedelta(int(delta_days))
        query = query.filter(Run.submitted >= start_date)
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
    site_id = request.form.get("site_id", None).upper()
    submission_id = request.form.get("submission_id")
    input_site_id = request.form.get("input_site_id", None).upper()
    input_endpoint = request.form.get("input_endpoint", None)
    output_endpoint = request.form.get("output_endpoint")
    output_dir = request.form.get("output_dir")
    wdl_file = request.form.get("wdl_file")
    json_file = request.form.get("json_file")
    tag = request.form.get("tag")
    compute_endpoint = config.conf.get_site(site_id, "globus_endpoint")
    globus = jaws_central.globus.GlobusService()

    if compute_endpoint is None:
        logger.error(
            f"Received run submission from {user} with invalid computing site ID: {site_id}"
        )
        abort(404, {"error": f'Unknown Site ID, "{site_id}"; try the "list-sites" command'})
    logger.info(f"User {user}: New run submission {submission_id} to {site_id}")

    # INSERT INTO RDB TO GET RUN ID
    run = Run(
        user_id=user,
        site_id=site_id,
        submission_id=submission_id,
        input_site_id=input_site_id,
        input_endpoint=input_endpoint,
        output_endpoint=output_endpoint,
        output_dir=output_dir,
        wdl_file=wdl_file,
        json_file=json_file,
        tag=tag,
        status="created",
    )
    try:
        db.session.add(run)
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Error inserting Run: {error}")
        abort(500, {"error": f"Error inserting Run into db: {error}"})
    try:
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        err_msg = f"Unable to insert new run in db: {error}"
        logger.exception(err_msg)
        abort(500, {"error": err_msg})
    logger.debug(f"User {user}: New run {run.id}")

    # Output directory is a subdirectory that includes the user id, site id and run id.
    # These are all placed in a common location with setgid sticky bits so that all
    # submitting users have access.
    output_dir += f"/{user}/{site_id}/{run.id}"
    src_host_path = config.conf.get_site(input_site_id, "globus_host_path")

    # We modify the output dir path since we know the endpoint of the returning source site. From here
    # a compute site can simply query the output directory and send.
    virtual_output_path = globus.virtual_transfer_path(output_dir, src_host_path)

    # Due to how the current database schema is setup, we have to update the output
    # directory from the model object itself immediately after insert.
    # TODO: Think of a better way to do this
    try:
        run.output_dir = virtual_output_path
    except Exception as error:
        db.session.rollback()
        err_msg = f"Unable to update output_dir in db: {error}"
        logger.exception(err_msg)
        abort(500, {"error": err_msg})
    try:
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        err_msg = f"Unable to update output_dir in db: {error}"
        logger.exception(err_msg)
        abort(500, {"error": err_msg})
    logger.debug(f"Updating output dir for run_id={run.id}")

    # SUBMIT GLOBUS TRANSFER
    manifest_file = request.files["manifest"]
    host_paths = {}
    host_paths["src"] = src_host_path
    host_paths["dest"] = config.conf.get_site(site_id, "globus_host_path")

    try:
        upload_task_id = globus.submit_transfer(
            f"Upload run {run.id}",
            host_paths,
            input_endpoint,
            compute_endpoint,
            manifest_file,
        )
    except globus_sdk.GlobusAPIError as error:
        run.status = "upload failed"
        db.session.commit()
        if error.code == "NoCredException":
            logger.warning(
                f"{user} submission {run.id} failed due to Globus {error.code}"
            )
            abort(
                401,
                {"error": error.message
                    + " -- Your access to the Globus endpoint has expired.  "
                    + "To reactivate, log-in to https://app.globus.org, go to Endpoints (left), "
                    + "search for the endpoint by name (if not shown), click on the endpoint, "
                    + "and use the button on the right to activate your credentials."},
            )
        else:
            logger.exception(
                f"{user} submission {run.id} failed for GlobusAPIError: {error}",
                exc_info=True,
            )
            abort(error.code, {"error": error.message})
    except globus_sdk.NetworkError as error:
        logger.exception(
            f"{user} submission {run.id} failed due to NetworkError: {error}",
            exc_info=True,
        )
        abort(500, {"error": f"Network Error: {error}"})
    except globus_sdk.GlobusError as error:
        logger.exception(
            f"{user} submission {run.id} failed for unknown error: {error}",
            exc_info=True,
        )
        abort(500, {"error": f"Unexpected error: {error}"})

    logger.debug(f"User {user}: Run {run.id} upload {upload_task_id}")

    # UPDATE RUN WITH UPLOAD TASK ID AND ADD LOG ENTRY
    run.upload_task_id = upload_task_id
    old_status = run.status
    new_status = run.status = "uploading"
    log = Run_Log(
        run_id=run.id,
        status_to=new_status,
        status_from=old_status,
        timestamp=run.updated,
        reason=f"upload_task_id={upload_task_id}",
    )
    try:
        db.session.add(log)
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Error inserting run log for Run {run.id}: {error}")
    try:
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        err_msg = f"Error inserting run log for Run {run.id}: {error}"
        logger.exception(err_msg)

    # GET CURRENT USER INFO
    try:
        current_user = db.session.query(User).get(user)
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, {"error": f"Db error; {e}"})

    # SEND TO SITE
    params = {
        "run_id": run.id,
        "user_id": user,
        "email": current_user.email,
        "submission_id": submission_id,
        "upload_task_id": upload_task_id,
        "output_endpoint": output_endpoint,
        "output_dir": output_dir,
    }
    a_site_rpc_client = rpc_index.rpc_index.get_client(run.site_id)
    logger.debug(f"User {user}: submit run: {params}")
    try:
        result = a_site_rpc_client.request("submit", params)
    except Exception as error:
        reason = f"RPC submit failed: {error}"
        logger.exception(reason)
        _submission_failed(user, run, reason)
        abort(500, {"error": reason})
    if "error" in result:
        reason = f"Error sending new run to {site_id}: {result['error']['message']}"
        logger.error(reason)
        _submission_failed(user, run, reason)
        abort(result["error"]["code"], {"error": result["error"]["message"]})

    # DONE
    result = {
        "run_id": run.id,
        "status": run.status,
        "site_id": site_id,
        "output_dir": output_dir,
        "tag": tag,
    }
    logger.info(f"User {user}: New run: {result}")
    return result, 201


def _submission_failed(user, run, reason):
    """Cancel upload and update run status"""
    _cancel_transfer(user, run.upload_task_id, run.id)
    _update_run_status(run, "submission failed", reason)


def _update_run_status(run, new_status, reason=None):
    """Update run table and insert run_logs entry."""
    status_from = run.status
    run.status = new_status
    try:
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Error updating run status in db: {error}")
    log = Run_Log(
        run_id=run.id,
        status_from=status_from,
        status_to=run.status,
        timestamp=run.updated,
        reason=reason,
    )
    try:
        db.session.add(log)
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Error insert run log entry: {error}")


def _get_run(user, run_id):
    """Return a Run object if found and accessible by current user, else abort.

    :param run_id: Run primary key
    :type run_id: str
    :return: the run ORM object for the record
    :rtype: sqlalchemy.model
    """
    try:
        run = db.session.query(Run).get(run_id)
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
            404, {"error": "No data available as the Run hasn't been submitted to Cromwell yet."}
        )


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
    run = _get_run(user, run_id)
    logger.info(f"User {user}: Get status of Run {run.id}")
    is_admin = _is_admin(user)
    info = _run_info(run, is_admin)
    return info, 200


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
    return _rpc_call(user, run_id, "get_task_status")


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
            db.session.query(Run_Log)
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
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "get_task_log")


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
    return _rpc_call(user, run_id, "run_metadata")


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
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "get_errors")


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

    # get run record
    run = _get_run(user, run_id)
    status = run.status

    # check if run can be cancelled
    if status == "cancelled":
        abort(400, {"error": "That Run had already been cancelled"})
    elif status == "download complete":
        abort(400, {"error": "It's too late to cancel; run is finished."})

    # mark as cancelled
    _cancel_run(run)

    # cancel active transfers
    if status.startswith("upload"):
        _cancel_transfer(user, run.upload_task_id, run_id)
    elif status.startswith("download"):
        _cancel_transfer(user, run.download_task_id, run_id)

    # tell Site to cancel
    return _rpc_call(user, run_id, "cancel_run")


def _cancel_run(run, reason="Cancelled by user"):
    """Update database record."""
    logger.debug(f"Run {run.id}: updating runs table to cancelled")
    status_from = run.status
    run.status = "cancelled"
    run.result = "cancelled"
    try:
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Error while updating run to 'cancelled': {error}")
    log = Run_Log(
        run_id=run.id,
        status_from=status_from,
        status_to=run.status,
        timestamp=run.updated,
        reason=reason,
    )
    try:
        db.session.add(log)
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.error(
            f"Error while adding run log entry to cancel run {run.id}: {error}"
        )


def _cancel_transfer(user: str, transfer_task_id: str, run_id: int) -> None:
    """Cancel a Globus transfer.

    :param user: user id
    :type user: str
    :param transfer_task_id: Globus transfer task id
    :type transfer_task_id: str
    :return: None; aborts on error.
    """
    logger.debug(f"Run {run_id}: Cancel transfer {transfer_task_id}")
    try:
        transfer_client = authorize_transfer_client()
    except globus_sdk.GlobusAPIError as error:
        logger.error(f"Error getting Globus transfer client: {error}")
        abort(500, {"error": "Globus error: {error}"})
    try:
        transfer_response = transfer_client.cancel_task(transfer_task_id)
        logger.debug(
            f"User {user} cancel upload {transfer_task_id} for run {run_id}: {transfer_response}"
        )
    except globus_sdk.GlobusAPIError as error:
        logger.exception(
            f"Failed to cancel {user}'s Globus transfer, {transfer_task_id}: {error}",
            exc_info=True,
        )
        abort(500, {"error": f"Globus error: {error}"})


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
        queue = (
            db.session.query(Run)
            .filter(and_(Run.user_id == user, Run.status.in_(run_active_states)))
            .all()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, {"error": f"Db error; {error}"})
    result = {}
    for run in queue:
        status = run.status
        if status == "cancelled":
            result[run.id] = "Already cancelled"
            continue
        elif status == "download complete":
            result[run.id] = "Too late to cancel; download complete"
            continue
        _cancel_run(run)
        if status.startswith("upload"):
            _cancel_transfer(user, run.upload_task_id, run.id)
        elif status.startswith("download"):
            _cancel_transfer(user, run.download_task_id, run.id)
        _rpc_call(user, run.id, "cancel_run")
        result[run.id] = "Cancelled"
    return result, 201


def authorize_transfer_client():
    """
    Create a globus transfer client using client id and client secret for credentials. More information
    can be found via Globus documentation:

    https://globus-sdk-python.readthedocs.io/en/stable/examples/client_credentials.html?highlight=secret

    :return: globus_sdk.TransferClient
    """
    client_id = config.conf.get("GLOBUS", "client_id")
    client_secret = config.conf.get("GLOBUS", "client_secret")
    try:
        client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
    except globus_sdk.GlobusAPIError as error:
        raise error
    scopes = "urn:globus:auth:scope:transfer.api.globus.org:all"
    try:
        authorizer = globus_sdk.ClientCredentialsAuthorizer(client, scopes)
    except globus_sdk.GlobusAPIError as error:
        raise error
    return globus_sdk.TransferClient(authorizer=authorizer)
