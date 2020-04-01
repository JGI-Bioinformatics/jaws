"""
Analysis (AKA Run) REST endpoints.
"""

import logging
from datetime import datetime, timedelta
from flask import abort, request
import globus_sdk
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config
from jaws_central.rpc_manager import rpc
from jaws_central.models import db, Run, User


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
        abort(404, "Run not found")
    if run.user_id != user:
        abort(401, "Access denied")
    site_rpc_call = rpc.get_client(run.site_id)
    params["cromwell_id"] = run.cromwell_id
    result = site_rpc_call.request(method, params)
    if "error" in result:
        abort(result["error"]["code"], result["error"]["message"])
    return result, 200


def user_queue(user):
    """Return the current user's unfinished runs.

    :param user: current user's ID
    :type user: str
    :return: details about any current runs
    :rtype: list
    """
    try:
        q = (
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
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, f"Db error; {e}")
    result = []
    for run in q:
        result.append(
            [run.id, run.submitted, run.status, run.submission_id, run.upload_task_id]
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
        q = (
            db.session.query(Run)
            .filter_by(user_id=user)
            .filter(Run.submitted >= start_date)
            .all()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, f"Db error; {e}")
    result = []
    for run in q:
        result.append(
            [run.id, run.submitted, run.status, run.submission_id, run.upload_task_id]
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
        abort(404, "Invalid Site ID")
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
        abort(404, f"Invalid Site ID: {site_id}")
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
        abort(500, "Unable to get Globus transfer client")
    tdata = globus_sdk.TransferData(
        transfer_client,
        input_endpoint,
        compute_endpoint,
        label=str(run.id),
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
    except globus_sdk.GlobusAPIError as e:
        run.status = "upload failed"
        db.session.commit()
        logger.exception(
            f"Error submitting Globus transfer for user {user}", exc_info=True
        )
        abort(500, f"Globus error: {e}")
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
        q = db.session.query(Run).get(run_id)
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, f"Db error; {e}")
    if not q:
        abort(404, "Run not found")
    return q


def _pre_cromwell(run):
    """Return True if the run hasn't been submitted to Cromwell yet.

    :param run: SQLAlchemy Model Run object
    :type param: sqlalchemy.model
    :return: True if pre-cromwell submission; false otherwise.
    :rtype: boolean
    """
    if run.status.startswith("upload") or run.status == "created":
        return True
    else:
        return False


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
    q = _get_run(run_id)
    result = {"status": q.status}
    return result, 200


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
    q = _get_run(run_id)
    if _pre_cromwell(q):
        return {"status": q.status}, 200
    else:
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
    q = _get_run(run_id)
    if _pre_cromwell(q):
        return {"status": q.status}, 200
    else:
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
    q = _get_run(run_id)
    if _pre_cromwell(q):
        return {"status": q.status}, 200
    else:
        return _rpc_call(user, run_id, "run_logs")


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
    q = _get_run(run_id)
    if _pre_cromwell(q):
        return {"status": q.status}, 200
    else:
        return _rpc_call(user, run_id, "failure_logs")


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
    q = _get_run(run_id)
    if _pre_cromwell(q):
        q.status = "canceled"
        db.session.commit()
        result = {"cancel": "OK"}
        return result, 200
    return _rpc_call(user, run_id, "cancel_run")
