import os
import logging
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from jaws_central import models
from jaws_central import config
from jaws_central import tasks
from datetime import datetime, timedelta
from jaws_central import jaws_constants
from sqlalchemy.exc import SQLAlchemyError
from jaws_rpc import rpc_index
from jaws_central.models_fsa import db, Run, User, Run_Log
from jaws_central.data_transfer import (
    DataTransferFactory,
    DataTransferProtocol,
    SiteTransfer,
    DataTransferAPIError,
    DataTransferNetworkError,
    DataTransferError,
)

logger = logging.getLogger(__package__)

#ECCE TODO

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
        run = db.session.query(Run).get(run_id)
    except SQLAlchemyError as e:
        logger.error(e)
        raise
    if not run:
        raise RunNotFoundError("Run not found; please check your run_id")
    if run.user_id != user and not _is_admin(user):
        raise RunAccessDeniedError(
            "Access denied; you cannot access to another user's workflow"
        )
    a_site_rpc_client = rpc_index.rpc_index.get_client(run.site_id)
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



#/ECCE


def file_size(file_path: str):
    """
    Checks if a file exists and is a file; if it does, check it's size.

    :param file_path: path to file
    :type file_path: str
    :return: Size in bytes if file exists; else None
    :rtype: int
    """
    if os.path.isfile(file_path):
        return os.path.getsize(file_path)
    else:
        return None


class RunDbError(Exception):
    pass


class RunNotFound(Exception):
    pass


class DataError(Exception):
    pass


class Run:
    """Class representing a single Run"""

    def __init__(self, session, **kwargs):
        self.session = session
        self.operations = {
            "uploading": self.check_if_upload_complete,
            "upload complete": self.submit_run,
            "succeeded": self.transfer_results,
            "failed": self.transfer_results,
            "downloading": self.check_if_download_complete,
        }

        if "submission_id" in kwargs:
            self._insert_run(**kwargs)
        elif "run_id" in kwargs:
            self._select_run(kwargs["run_id"])
        elif "model" in kwargs:
            model = kwargs["model"]
            assert model.id is not None
            self.model = model

# ECCE TODO
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



    def _insert_run(self, **kwargs) -> None:
        """Insert run record into rdb"""
        try:
            self.model = models.Run(
                id=int(kwargs["run_id"]),
                user_id=kwargs["user_id"],
                email=kwargs["email"],
                submission_id=kwargs["submission_id"],
                upload_task_id=kwargs["upload_task_id"],
                output_endpoint=kwargs["output_endpoint"],
                output_dir=kwargs["output_dir"],
                status="uploading",
            )
        except SQLAlchemyError as error:
            raise (f"Error creating model for new Run {kwargs['run_id']}: {error}")
        try:
            self.session.add(self.model)
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            raise (error)

    def _select_run(self, run_id):
        """Select run record from rdb"""
        try:
            self.model = self.session.query(models.Run).get(run_id)
        except IntegrityError as error:
            logger.warn(f"Run {run_id} not found: {error}")
            raise RunNotFound(f"Run {run_id} not found")
        except SQLAlchemyError as error:
            err_msg = f"Unable to select run, {run_id}: {error}"
            logger.error(err_msg)
            raise RunDbError(err_msg)

    @property
    def status(self) -> str:
        """Return the current state of the run."""
        return self.model.status

    def check_status(self) -> None:
        """Check the run's status, promote to next state if ready"""
        status = self.model.status
        if status in self.operations:
            return self.operations[status]()

    def cancel(self) -> None:
        """Cancel a run, instruct Cromwell to abort if required."""
        self.update_run_status("cancelled")

        status = self.model.status
        if status == "cancelled":
            abort(400, {"error": "That Run had already been cancelled"})
        elif status == "download complete":
            abort(400, {"error": "It's too late to cancel; run is finished."})
        if self.model.cromwell_run_id and self.model.status in [
            "submitted",
            "queued",
            "running",
        ]:
            # the Run is in a cancellable state
            try:
                rpc_cancel(self.model.id)  # TODO
            except Exception as error:
                logger.warn(f"Unable to cancel Run {self.model.id}: {error}")
                raise
        # TODO
        cancelled = _cancel_run(user, run)
        return {run_id: cancelled}, 201


    # TODO
    def _cancel_run(self, user, run, reason="Cancelled by user"):
        """
        Cancel a Run.

        :param run: Run SqlAlchemy ORM object
        :type run: obj
        """
        transfer_type = SiteTransfer.type[run.site_id.upper()]
        data_transfer = DataTransferFactory(transfer_type)
        status = run.status

        if status.startswith("upload"):
            data_transfer.cancel_transfer(run.upload_task_id)
        elif status.startswith("download"):
            data_transfer.cancel_transfer(run.upload_task_id)
        try:
            _rpc_call(user, run.id, "cancel_run")
        except Exception as error:
            logger.error(f"Error canceling run {run.id}: {error}")
            # ignore error, cancel anyway
        try:
            _update_run_status(run, "cancelled", reason)
        except Exception as error:
            return f"cancel failed; {error}"
        else:
            return "cancelled"









    def metadata(self) -> str:
        """
        Get metadata from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.model.cromwell_run_id:
            return cromwell.get_metadata(self.model.cromwell_run_id).data
        else:
            return None

    def outputs(self) -> str:
        """
        Get outputs from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.model.cromwell_run_id:
            metadata = cromwell.get_metadata(self.model.cromwell_run_id)
            return metadata.outputs(relpath=True)
        else:
            return None

    def errors(self) -> str:
        """Get errors report from Cromwell service"""
        if self.model.cromwell_run_id:
            metadata = cromwell.get_metadata(self.model.cromwell_run_id)
            return metadata.errors()
        else:
            return None

    def upload_status(self, data_transfer: DataTransferProtocol) -> str:
        try:
            status = data_transfer.transfer_status(self.model.upload_task_id)
        except DataTransferError as error:
            logger.exception(
                f"Failed to check upload {self.model.upload_task_id}: {error}"
            )
        else:
            return status

    def download_status(self, data_transfer: DataTransferProtocol) -> str:
        try:
            status = data_transfer.transfer_status(self.model.download_task_id)
        except DataTransferError as error:
            logger.exception(
                f"Failed to check download {self.model.download_task_id}: {error}"
            )
        else:
            return status

    def check_if_upload_complete(self) -> None:
        """
        If the upload is complete, update the run's status.
        """
        logger.debug(f"Run {self.model.id}: Check upload status")
        data_transfer_type = self._get_data_transfer_type()
        data_transfer = DataTransferFactory(data_transfer_type)
        status = self.upload_status(data_transfer)
        if status == SiteTransfer.status.failed:
            self.update_run_status("upload failed")
        elif status == SiteTransfer.status.inactive:
            self.update_run_status(
                "upload inactive", "The JAWS endpoint authorization has expired"
            )
        elif status == SiteTransfer.status.succeeded:
            self.update_run_status("upload complete")

    def uploads_file_path(self):
        uploads_dir = config.conf.get("SITE", "uploads_dir")
        return os.path.join(uploads_dir, self.model.user_id, self.model.submission_id)

    def get_run_input(self) -> list:
        """
        Check if the input files are valid and return a list of their file handles.
        The list contains wdl_file, json_file, and optionally zip_file.
        Raise on error.

        :return: list of [wdl,json,zip] file paths
        :rtype: list
        """
        infile_parameters = [
            ("wdl", True, False),
            ("json", True, False),
            ("zip", False, True),
            ("options.json", False, False),
        ]
        file_handles = []
        file_path = self.uploads_file_path()
        for (suffix, required, is_binary) in infile_parameters:
            a_file = f"{file_path}.{suffix}"
            a_file_size = file_size(a_file)
            if required:
                if a_file_size is None:
                    raise DataError(f"Input {suffix} file not found: {a_file}")
                elif a_file_size == 0:
                    raise DataError(f"Input {suffix} file is 0-bytes: {a_file}")
            elif a_file_size is not None and a_file_size == 0:
                raise DataError(f"Input {suffix} file is 0-bytes: {a_file}")
            a_fh = None
            filetype = "rb" if is_binary else "r"
            if a_file_size:
                try:
                    a_fh = open(a_file, filetype)
                except IOError:
                    raise
            file_handles.append(a_fh)
        return file_handles

    def submit_run(self) -> None:
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {self.model.id}: Submit to Site")
        # NYI
        #self.update_run_status("sent")

# ECCE TODO
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

    if compute_endpoint is None:
        logger.error(
            f"Received run submission from {user} with invalid computing site ID: {site_id}"
        )
        abort(
            404,
            {"error": f'Unknown Site ID, "{site_id}"; try the "list-sites" command'},
        )
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

    # Setup data transfer object. The data transfer using either globus or AWS is based on site_id. This is defined
    # in the SiteTransfer.type var.
    data_transfer_type = SiteTransfer.type[site_id.upper()]
    data_transfer = DataTransferFactory(data_transfer_type)
    metadata = {"label": f"Upload run {run.id}"}

    if data_transfer_type == "globus_transfer":
        src_host_path = config.conf.get_site(input_site_id, "globus_host_path")
        metadata["host_paths"] = {
            "src": src_host_path,
            "dest": config.conf.get_site(site_id, "globus_host_path"),
        }
        metadata["input_endpoint"] = input_endpoint
        metadata["compute_endpoint"] = compute_endpoint
        metadata["run_id"] = run.id

        # We modify the output dir path since we know the endpoint of the returning source site. From here
        # a compute site can simply query the output directory and send.
        virtual_output_path = data_transfer.virtual_transfer_path(
            output_dir, src_host_path
        )
    else:
        virtual_output_path = (
            output_dir  # not sure what the virtual path for AWS should be ???
        )

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

    # SUBMIT FILE TRANSFER
    manifest_files = request.files["manifest"]

    logger.debug(f"Transferring files using {data_transfer_type}")

    try:
        upload_task_id = data_transfer.submit_transfer(metadata, manifest_files)
    except DataTransferAPIError as error:
        run.status = "upload failed"
        db.session.commit()
        if error.code == "NoCredException":
            logger.warning(
                f"{user} submission {run.id} failed due to Globus {error.code}"
            )
            abort(
                401,
                {
                    "error": error.message
                    + " -- Your access to the Globus endpoint has expired.  "
                    + "To reactivate, log-in to https://app.globus.org, go to Endpoints (left), "
                    + "search for the endpoint by name (if not shown), click on the endpoint, "
                    + "and use the button on the right to activate your credentials."
                },
            )
        else:
            logger.exception(
                f"{user} submission {run.id} failed for GlobusAPIError: {error}",
                exc_info=True,
            )
            abort(error.code, {"error": error.message})
    except DataTransferNetworkError as error:
        logger.exception(
            f"{user} submission {run.id} failed due to NetworkError: {error}",
            exc_info=True,
        )
        abort(500, {"error": f"Network Error: {error}"})
    except DataTransferError as error:
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
        _submission_failed(user, run, reason, data_transfer)
        abort(500, {"error": reason})
    if "error" in result:
        reason = f"Error sending new run to {site_id}: {result['error']['message']}"
        logger.error(reason)
        _submission_failed(user, run, reason, data_transfer)
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

def _submission_failed(
    user: str, run: str, reason: str, data_transfer: DataTransferProtocol
):
    """Cancel upload and update run status"""
    data_transfer.cancel_transfer(run.upload_task_id)
    _update_run_status(run, "submission failed", reason)


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

#/ECCE


    def check_run_cromwell_status(self) -> None:
        """
        Check Cromwell for the status of the Run.
        """
        logger.debug(f"Run {self.model.id}: Check Cromwell status")
        # TODO RPC

    def _get_data_transfer_type(self) -> str:
        """Lookup the site id to determine what type of data transfer method needs to be used (e.g., globus, aws).

        :return: data transfer type based on site id (e.gl, 'globus_transfer', 'aws_transfer')
        :rtype: str
        """
        site_id = config.conf.get("SITE", "id")
        transfer_type = SiteTransfer.type[site_id.upper()]
        return transfer_type

    def transfer_results(self) -> None:
        """
        Send run output via data transfer
        """
        logger.debug(f"Run {self.model.id}: Download output")
        try:
            metadata = cromwell.get_metadata(self.model.cromwell_run_id)
        except CromwellServiceError as error:
            # if Cromwell service not available, do not fail run and try again later
            logger.debug(f"Cromwell service error: {error}")
            return
        except CromwellRunError as error:
            # if there's something wrong with this particular (failed) run, promote to next state
            logger.debug(f"Run {self.model.id}: {error}")
            self.update_run_status("download complete", f"Cromwell run error: {error}")
            return
        cromwell_workflow_dir = metadata.workflow_root()
        if not cromwell_workflow_dir:
            # This run failed before a folder was created; nothing to xfer
            self.update_run_status("download complete", "No run folder was created")
            return

        data_transfer_type = self._get_data_transfer_type()
        data_transfer = DataTransferFactory(data_transfer_type)
        metadata = {
            "label": f"Run {self.model.id}"
        }

        if data_transfer_type == 'globus_transfer':
            metadata['dest_endpoint'] = self.model.output_endpoint

        logger.debug(f"Transferring files using {data_transfer_type}")

        try:
            transfer_task_id = self.transfer_files(data_transfer, metadata, cromwell_workflow_dir,
                                                   self.model.output_dir)
        except DataTransferError as error:
            logger.error(f"error while submitting transfer: {error}")
            raise
        else:
            logger.debug(f"Transfer task id={transfer_task_id}")
            self.model.download_task_id = transfer_task_id
            self.update_run_status(
                "downloading", f"download_task_id={self.model.download_task_id}"
            )

    def transfer_files(self, data_transfer: DataTransferProtocol, metadata: dict, src_dir: str, dst_dir: str) -> str:
        """Transfer files using the data_transfer object (e.g., via globus or aws)."""
        return data_transfer.submit_download(metadata, src_dir, dst_dir)

    def check_if_download_complete(self) -> None:
        """
        If download is complete, change state.
        """
        logger.debug(f"Run {self.model.id}: Check download status")
        data_transfer_type = self._get_data_transfer_type()
        data_transfer = DataTransferFactory(data_transfer_type)
        status = self.download_status(data_transfer)
        if status == SiteTransfer.status.succeeded:
            self.update_run_status("download complete")
        elif status == SiteTransfer.status.failed:
            self.update_run_status("download failed")

    def update_run_status(self, status_to, reason=None) -> None:
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        status_from = self.model.status
        logger.info(f"Run {self.model.id}: now {status_to}")
        timestamp = datetime.utcnow()
        self._update_run_status(status_to, timestamp)
        self._insert_run_log(status_from, status_to, timestamp, reason)

# ECCE TODO
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



    def _update_run_status(self, new_status, timestamp) -> None:
        """
        Update Run's current status in 'runs' table
        """
        try:
            self.model.status = new_status
            self.model.updated = timestamp
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Unable to update Run {self.model.id}: {error}")
            raise

    def _insert_run_log(self, status_from, status_to, timestamp, reason) -> None:
        """
        Save record of state transition in 'run_logs' table
        """
        try:
            log_entry = models.Run_Log(
                run_id=self.model.id,
                status_from=status_from,
                status_to=status_to,
                timestamp=timestamp,
                reason=reason,
            )
            self.session.add(log_entry)
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(
                f"Failed to insert log for Run {self.model.id} ({status_to}): {error}"
            )
            raise

    # ECCE TODO
    def task_status(self):
        _abort_if_pre_cromwell(run)
        return rpc_call(user, run_id, "get_task_status")


    def run_log(self):
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


    def task_log(self):
        _abort_if_pre_cromwell(run)
        return rpc_call(user, run_id, "get_task_log")


    def task_summary(self):
        _abort_if_pre_cromwell(run)
        return rpc_call(user, run_id, "get_task_summary")


    def metadata(self):
        _abort_if_pre_cromwell(run)
        return rpc_call(user, run_id, "run_metadata")


    def run_outputs(self):
        _abort_if_pre_cromwell(run)
        return rpc_call(user, run_id, "run_outputs")


    def errors(self) -> dict:
        _abort_if_pre_cromwell(run)
        return rpc_call(user, run_id, "get_errors")

    # /ECCE


def check_active_runs(session) -> None:
    """
    Get active runs from db and have each check and update their status.
    """
    rows = _select_active_runs(session)
    for model in rows:
        run = Run(session, model=model)
        run.check_status()


def _select_active_runs(session) -> list:
    """
    Select runs in particular states from db.
    """
    active_states = [
        "uploading",
        "upload complete",
        "submitted",
        "queued",
        "running",
        "succeeded",
        "failed",
        "downloading",
    ]
    try:
        rows = (
            session.query(models.Run).filter(models.Run.status.in_(active_states)).all()
        )
    except SQLAlchemyError as error:
        logger.warning(f"Failed to select active runs from db: {error}", exc_info=True)
        return []
    return rows


# ECCE TODO
def _run_info(run, is_admin: bool = False, verbose: bool = False):
    """
    Given a SQLAlchemy model for a Run, create a dict with the desired fields.
    :param run: Run object
    :type run: model
    :param is_admin: True if current user is an administrator
    :type is_admin: bool
    :param verbose: True if all fields desired
    :type verbose: bool
    :return: selected fields
    :rtype: dict
    """
    info = {}
    complete = True if (is_admin or verbose) else False
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
    for run in rows:
        runs.append(_run_info(run, is_admin))
    return runs, 200


def _select_runs(user: str, **kwargs):
    """Select runs from db.

    :param user: current user's ID
    :type user: str
    :return: Runs matching search criteria
    :rtype: list
    """
    query = db.session.query(Run)
    if "all_users" in kwargs and kwargs["all_users"] is True:
        pass
    else:
        query = query.filter(Run.user_id == user)
    if "active_only" in kwargs and kwargs["active_only"] is True:
        query = query.filter(Run.status.in_(run_active_states))
    if "site_id" in kwargs:
        site_id = kwargs["site_id"].upper()
        if site_id != "ALL":
            query = query.filter(Run.site_id == site_id)
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


# /ECCE
