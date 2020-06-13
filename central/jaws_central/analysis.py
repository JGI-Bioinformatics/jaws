"""
Analysis (AKA Run) REST endpoints.
"""

import logging
from datetime import datetime, timedelta
from flask import abort, request
import globus_sdk
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config
from jaws_central import jaws_constants
from jaws_rpc import rpc_index
from jaws_central.models import db, Run, User, Run_Log


logger = logging.getLogger(__package__)


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
        abort(401, "Access denied; you cannot access to another user's workflow")
    site_rpc_call = rpc_index.index.get_client(run.site_id)
    params["user"] = user
    params["cromwell_id"] = run.cromwell_id
    logger.info(f"User {user} RPC {method} params {params}")
    result = site_rpc_call.request(method, params)
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
    try:
        queue = (
            db.session.query(Run)
            .filter_by(user_id=user)
            .filter(
                Run.status.in_(
                    [
                        "created",
                        "uploading",
                        "queued",
                        "running",
                        "post-processing",
                        "downloading",
                    ]
                )
            )
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
                "cromwell_id": run.cromwell_id,
                "status": run.status,
                "site_id": run.site_id,
                "submitted": run.submitted,
                "updated": run.updated,
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
    try:
        history = (
            db.session.query(Run)
            .filter_by(user_id=user)
            .filter(Run.submitted >= start_date)
            .all()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, f"Db error; {e}")
    result = []
    for run in history:
        result.append(
            {
                "id": run.id,
                "submission_id": run.submission_id,
                "cromwell_id": run.cromwell_id,
                "status": run.status,
                "site_id": run.site_id,
                "submitted": run.submitted,
                "updated": run.updated,
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
        abort(404, f"Unknown Site ID; {site_id} is not one of our sites")
    logger.info(f"New submission from {user} to {site_id}")

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
    db.session.add(run)
    db.session.commit()

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
            500,
            "User Globus error; have you granted JAWS access via the 'login' command?",
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
    logger.info(f"Upload for {user} from {input_endpoint} to {compute_endpoint}")
    manifest_file = request.files["manifest"]
    manifest = manifest_file.read().splitlines()
    for line in manifest:
        line = line.decode("UTF-8")
        [source_path, dest_path, inode_type] = line.split("\t")
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
        logger.info(f"{user} submission {run.id} failed due to Globus {error.code}")
        if error.code == "NoCredException":
            abort(
                error.http_status,
                error.message
                + " -- Your access to the Globus endpoint has expired.  "
                + "To reactivate, log-in to https://globus.org, go to Endpoints (left), "
                + "search for the endpoint by name (if not shown), click on the endpoint, "
                + "and use the button on the right to activate your credentials.",
            )
        else:
            logger.exception(
                f"{user} submission {run.id} failed for GlobusAPIError", exc_info=True
            )
            abort(error.http_status, error.message)
    except globus_sdk.NetworkError as error:
        logger.info(
            f"{user} submission {run.id} failed due to NetworkError", exc_info=True
        )
        abort(503, f"Network Error: {error}")
    except globus_sdk.GlobusError as error:
        logger.exception(
            f"{user} submission {run.id} failed for unknown error", exc_info=True
        )
        abort(500, f"Unexpected error: {error}")
    upload_task_id = transfer_result["task_id"]

    # UPDATE RUN WITH UPLOAD TASK ID AND RETURN RESULTS
    run.upload_task_id = upload_task_id
    run.status = "uploading"
    db.session.commit()
    result = {
        "run_id": run.id,
        "submission_id": submission_id,
        "status": run.status,
        "upload_task_id": upload_task_id,
        "site_id": site_id,
        "output_endpoint": output_endpoint,
        "output_dir": output_dir,
    }
    return result, 201


def _get_run(run_id):
    """Return a Run object if found, else abort.

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
    return run


def __abort_if_pre_cromwell(run):
    """Returns if run was submitted to Cromwell, aborts otherwise.

    :param run: SQLAlchemy Model Run object
    :type param: sqlalchemy.model
    :return: True if pre-cromwell submission; false otherwise.
    :rtype: boolean
    """
    if run.status.startswith("upload") or run.status == "created":
        abort(204, "No data; this run hasn't begun execution yet.")


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
    run = _get_run(run_id)
    result = {
        "id": run.id,
        "submission_id": run.submission_id,
        "cromwell_id": run.cromwell_id,
        "status": run.status,
        "status_detail": jaws_constants.run_status_msg.get(run.status, ""),
        "site_id": run.site_id,
        "submitted": run.submitted,
        "updated": run.updated,
        "input_site_id": run.input_site_id,
        "input_endpoint": run.input_endpoint,
        "upload_task_id": run.upload_task_id,
        "output_endpoint": run.output_endpoint,
        "output_dir": run.output_dir,
        "download_task_id": run.download_task_id,
    }
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
    try:
        query = (
            db.session.query(Run_Log)
            .filter_by(run_id=run_id)
            .order_by(Run_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.error(error)
        abort(500, f"Db error; {error}")
    table = []
    for log in query:
        row = [
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            log.reason,
        ]
        table.append(row)
    return table, 200


def task_status(user, run_id):
    """
    Retrieve run status with task-level detail.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: The status of each task in a run, if found; abort otherwise
    :rtype: dict
    """
    run = _get_run(run_id)
    __abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "task_status")


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
    run = _get_run(run_id)
    __abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "run_metadata")


def run_logs(user, run_id):
    """
    Retrieve the logs of a run.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell metadata for the run, if found; abort otherwise
    :rtype: dict
    """
    run = _get_run(run_id)
    __abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "run_logs")


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
    run = _get_run(run_id)
    __abort_if_pre_cromwell(run)
    return _rpc_call(user, run_id, "output", {"failed_only": False})


def failure_logs(user, run_id):
    """
    Retrieve the logs for failed tasks.

    :param user: current user's ID
    :type user: str
    :param run_id: unique identifier for a run
    :type run_id: int
    :return: Cromwell stdout/stderr/submit files for failed tasks.
    :rtype: str
    """
    run = _get_run(run_id)
    __abort_if_pre_cromwell(run)
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
    run = _get_run(run_id)
    if run.user_id != user:
        abort(401, "Access denied; you cannot cancel another user's runs")
    # the run state will be changed to "canceled" regardless of result of any other actions below
    status = run.status
    run.status = "canceled"
    db.session.commit()
    logger.info(f"{user} cancelled run {run_id} in {run.status} state")
    if status == "submitted" or status == "running":
        # if running, instruct Cromwell to cancel
        params = {"cromwell_id": run.cromwell_id, "user": user}
        _rpc_call(user, run_id, "cancel_run", params)  # will abort if fail
    elif status == "uploading":
        _cancel_transfer(user, run.upload_task_id, run_id)  # will abort if fail
    elif status == "downloading":
        _cancel_transfer(user, run.download_task_id, run_id)
    result = {"cancel": "OK"}
    return result, 201


def _cancel_transfer(user: str, transfer_task_id: str, run_id: int) -> None:
    """Cancel a Globus transfer.

    :param user: user id
    :type user: str
    :param transfer_task_id: Globus transfer task id
    :type transfer_task_id: str
    :return: None; aborts on error.
    """
    logger.info(f"{user} cancelled Globus transfer {transfer_task_id} for run {run_id}")
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
        logger.exception(
            f"Error getting transfer client for user {user}", exc_info=True
        )
        abort(500, "Unable to get Globus transfer client; perhaps try 'login' command")
    try:
        transfer_response = transfer_client.cancel_task(transfer_task_id)
        logger.info(
            f"User {user} cancel upload {transfer_task_id} for run {run_id}: {transfer_response}"
        )
    except globus_sdk.GlobusAPIError as e:
        logger.exception(
            f"Failed to cancel {user}'s Globus transfer, {transfer_task_id}: {transfer_response}",
            exc_info=True,
        )
        abort(500, f"Globus error: {e}")
