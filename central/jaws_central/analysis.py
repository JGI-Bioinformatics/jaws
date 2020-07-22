"""
Analysis (AKA Run) REST endpoints.
"""

import logging
from datetime import datetime, timedelta
from flask import abort, request
import globus_sdk
from sqlalchemy.exc import SQLAlchemyError
import collections
from jaws_central import config
from jaws_central import jaws_constants
from jaws_rpc import rpc_index
from jaws_central.models_fsa import db, Run, User, Run_Log, Job_Log


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
    "downloading"
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
        abort(500, f"Db error; {e}")
    if not run:
        abort(404, "Run not found; please check your run_id")
    if run.user_id != user:
        try:
            current_user = db.session.query(User).get(user)
        except SQLAlchemyError as error:
            logger.error(error)
            abort(500, f"Db error; {error}")
        if not current_user.is_admin:
            abort(401, "Access denied; you cannot access to another user's workflow")
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
        abort(result["error"]["code"], result["error"]["message"])
    return result["result"], 200


def user_queue(user):
    """Return the current user's unfinished runs.

    :param user: current user's ID
    :type user: str
    :return: details about any current runs
    :rtype: list
    """
    logger.info(f"User {user}: Get queue")
    try:
        queue = (
            db.session.query(Run)
            .filter_by(user_id=user)
            .filter(Run.status.in_(run_active_states))
            .all()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    result = []
    for run in queue:
        result.append(
            {
                "id": run.id,
                "submission_id": run.submission_id,
                "cromwell_run_id": run.cromwell_run_id,
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
            }
        )
    return result, 200


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
        history = (
            db.session.query(Run)
            .filter_by(user_id=user)
            .filter(Run.submitted >= start_date)
            .all()
        )
    except SQLAlchemyError as error:
        logger.exception(f"Failed to select run history: {error}")
        abort(500, f"Db error; {error}")
    result = []
    for run in history:
        result.append(
            {
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
            }
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
        result.append(site_id)
    return result, 200


def get_site(user, site_id):
    """Get parameters of a Site, required to submit a run.

    :param user: current user's ID
    :type user: str
    :param site_id: a JAWS-Site's ID
    :type site_id: str
    :return: globus endpoint id and staging path
    :rtype: dict
    """
    logger.debug(f"User {user}: Get info for site {site_id}")
    result = config.conf.get_site_info(site_id)
    if result is None:
        abort(404, f'Unknown Site ID; "{site_id}" is not one of our sites')
    result["staging_subdir"] = f'{result["staging_subdir"]}/{user}'
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
    compute_endpoint = config.conf.get_site(site_id, "globus_endpoint")
    if compute_endpoint is None:
        logger.error(
            f"Received run submission from {user} with invalid computing site ID: {site_id}"
        )
        abort(404, f'Unknown Site ID, "{site_id}"; try the "list-sites" command')
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
        status="created",
    )
    try:
        db.session.add(run)
    except Exception as error:
        logger.exception(f"Error inserting Run: {error}")
        abort(500, f"Error inserting Run into db: {error}")
    db.session.commit()
    logger.debug(f"User {user}: New run {run.id}")

    # SUBMIT GLOBUS TRANSFER
    try:
        current_user = db.session.query(User).get(user)
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, f"Db error; {e}")
    transfer_rt = current_user.transfer_refresh_token
    client_id = config.conf.get("GLOBUS", "client_id")
    try:
        client = globus_sdk.NativeAppAuthClient(client_id)
        authorizer = globus_sdk.RefreshTokenAuthorizer(transfer_rt, client)
        transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
    except globus_sdk.GlobusAPIError:
        run.status = "upload failed"
        db.session.commit()
        logger.exception(
            f"Error getting transfer client for user {user}", exc_info=True
        )
        abort(
            401,
            "Globus access denied; have you granted JAWS access via the 'login' command?",
        )
    tdata = globus_sdk.TransferData(
        transfer_client,
        input_endpoint,
        compute_endpoint,
        label=f"Upload run {run.id}",
        sync_level="checksum",
        verify_checksum=True,
        preserve_timestamp=True,
        notify_on_succeeded=False,
        notify_on_failed=True,
        notify_on_inactive=True,
        skip_activation_check=False,
    )
    manifest_file = request.files["manifest"]
    manifest = manifest_file.read().splitlines()
    for line in manifest:
        line = line.decode("UTF-8")
        source_path, dest_path, inode_type = line.split("\t")
        logger.debug(f"add transfer: {source_path} -> {dest_path}")
        if inode_type == "D":
            tdata.add_item(source_path, dest_path, recursive=True)
        else:
            tdata.add_item(source_path, dest_path, recursive=False)
    try:
        transfer_result = transfer_client.submit_transfer(tdata)
    except globus_sdk.GlobusAPIError as error:
        run.status = "upload failed"
        db.session.commit()
        if error.code == "NoCredException":
            logger.warning(f"{user} submission {run.id} failed due to Globus {error.code}")
            abort(
                401,
                error.message
                + " -- Your access to the Globus endpoint has expired.  "
                + "To reactivate, log-in to https://app.globus.org, go to Endpoints (left), "
                + "search for the endpoint by name (if not shown), click on the endpoint, "
                + "and use the button on the right to activate your credentials.",
            )
        else:
            logger.exception(
                f"{user} submission {run.id} failed for GlobusAPIError: {error}", exc_info=True
            )
            abort(error.code, error.message)
    except globus_sdk.NetworkError as error:
        logger.exception(
            f"{user} submission {run.id} failed due to NetworkError: {error}", exc_info=True
        )
        abort(500, f"Network Error: {error}")
    except globus_sdk.GlobusError as error:
        logger.exception(
            f"{user} submission {run.id} failed for unknown error: {error}", exc_info=True
        )
        abort(500, f"Unexpected error: {error}")
    upload_task_id = transfer_result["task_id"]
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
        logger.exception(f"Error inserting run log for Run {run.id}: {error}")
    db.session.commit()

    # SEND TO SITE
    params = {
        "run_id": run.id,
        "user_id": user,
        "email": current_user.email,
        "transfer_refresh_token": current_user.transfer_refresh_token,
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
        abort(500, reason)
    if "error" in result:
        reason = f"Error sending new run to {site_id}: {result['error']['message']}"
        logger.error(reason)
        _submission_failed(user, run, reason)
        abort(result["error"]["code"], result["error"]["message"])

    # DONE
    result = {
        "run_id": run.id,
        "submission_id": submission_id,
        "status": run.status,
        "upload_task_id": upload_task_id,
        "site_id": site_id,
        "output_endpoint": output_endpoint,
        "output_dir": output_dir,
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
    db.session.commit()
    log = Run_Log(
        run_id=run.id,
        status_from=status_from,
        status_to=run.status,
        timestamp=run.updated,
        reason=reason,
    )
    db.session.add(log)
    db.session.commit()


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
        abort(500, f"Db error; {error}")
    if not run:
        abort(404, "Run not found; please check your run_id")
    if run.user_id != user:
        try:
            current_user = db.session.query(User).get(user)
        except SQLAlchemyError as error:
            logger.error(error)
            abort(500, f"Db error; {error}")
        if not current_user.is_admin:
            abort(401, "Access denied; you are not the owner of that Run.")
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
            404, "No data available as the Run hasn't been submitted to Cromwell yet."
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
    result = {
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
    }
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
    run = _get_run(user, run_id)
    logger.info(f"User {user}: Get task-status of Run {run.id}")
    if not run:
        abort(404, f"Run {run_id} does not exist")

    # get job log entries, sorted by timestamp,
    # save in dict in which only the latest entry is kept
    try:
        query = (
            db.session.query(Job_Log)
            .filter_by(run_id=run_id)
            .order_by(Job_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting from job_log: {error}")
        abort(500, f"Db error; {error}")
    tasks = collections.OrderedDict()
    for log in query:
        task_name = log.task_name
        reason = log.reason if log.reason else ""
        tasks[task_name] = [
            log.task_name,
            log.attempt,
            log.cromwell_job_id,
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            reason,
            jaws_constants.task_status_msg.get(log.status_to, "")
        ]
        tasks.move_to_end(task_name)
    return list(tasks.values()), 200


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
        abort(500, f"Db error; {error}")
    table = []
    for log in query:
        reason = log.reason if log.reason else ""
        row = [
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            reason
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
    run = _get_run(user, run_id)
    logger.info(f"User {user}: Get task-log of Run {run.id}")
    try:
        query = (
            db.session.query(Job_Log)
            .filter_by(run_id=run_id)
            .order_by(Job_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting from job_logs: {error}")
        abort(500, f"Db error; {error}")
    table = []
    for log in query:
        reason = log.reason if log.reason else ""
        row = [
            log.task_name,
            log.attempt,
            log.cromwell_job_id,
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            reason
        ]
        table.append(row)
    return table, 200


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


def output(user, run_id):
    """
    Retrieve the stdout/stderr output of all Tasks.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: stdout/stderr file contents
    :rtype: str
    """
    logger.info(f"User {user}: Get output for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "output", {"failed_only": False})


def failed_output(user, run_id):
    """
    Retrieve the logs for failed tasks.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell stdout/stderr/submit files for failed tasks.
    :rtype: str
    """
    logger.info(f"User {user}: Get failed-task output for Run {run_id}")
    run = _get_run(user, run_id)
    _abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "output", {"failed_only": True})


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
        abort(400, "That Run had already been cancelled")
    elif status == "download complete":
        abort(400, "It's too late to cancel; run is finished.")

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
        logger.error(f"Error while adding run log entry to cancel run {run.id}: {error}")


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
        current_user = db.session.query(User).get(user)
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, f"Db error; {e}")
    transfer_rt = current_user.transfer_refresh_token
    client_id = config.conf.get("GLOBUS", "client_id")
    try:
        client = globus_sdk.NativeAppAuthClient(client_id)
        authorizer = globus_sdk.RefreshTokenAuthorizer(transfer_rt, client)
        transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
    except globus_sdk.GlobusAPIError as error:
        logger.exception(
            f"Error getting transfer client for user {user}: {error}", exc_info=True
        )
        abort(500, "Unable to get Globus transfer client; perhaps try 'login' command")
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
        abort(500, f"Globus error: {error}")
