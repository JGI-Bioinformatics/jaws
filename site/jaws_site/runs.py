"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import os
import logging
from datetime import datetime
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
import json
from jaws_site import models
from jaws_site import config
from jaws_site import tasks
from jaws_site.datatransfer_protocol import (
    DataTransferError,
    SiteTransfer,
    DataTransferProtocol,
    DataTransferFactory,
)
from jaws_site.cromwell import (
    Cromwell,
    CromwellError,
    CromwellServiceError,
    CromwellRunError,
)

logger = logging.getLogger(__package__)

cromwell = Cromwell(config.conf.get("CROMWELL", "url"))


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
            "submitted": self.check_run_cromwell_status,
            "queued": self.check_run_cromwell_status,
            "running": self.check_run_cromwell_status,
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
        if self.model.cromwell_run_id and self.model.status in [
            "submitted",
            "queued",
            "running",
        ]:
            try:
                cromwell.abort(self.model.cromwell_run_id)
            except CromwellError as error:
                logger.warn(f"Cromwell error cancelling Run {self.model.id}: {error}")
                raise

    def metadata(self) -> str:
        """
        Get metadata from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.model.cromwell_run_id:
            return cromwell.get_metadata(self.model.cromwell_run_id).data
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
        Check if the input files are valid and return a list of their paths.
        The list contains wdl_file, json_file, and optionally zip_file.
        Raise on error.

        :return: list of [wdl,json,zip] file paths
        :rtype: list
        """
        suffixes_and_required = [
            ("wdl", True),
            ("json", True),
            ("zip", False),
            ("options.json", False),
        ]
        files = []
        file_path = self.uploads_file_path()
        for (suffix, required) in suffixes_and_required:
            a_file = f"{file_path}.{suffix}"
            a_file_size = file_size(a_file)
            if required:
                if a_file_size is None:
                    raise DataError(f"Input {suffix} file not found: {a_file}")
                elif a_file_size == 0:
                    raise DataError(f"Input {suffix} file is 0-bytes: {a_file}")
            elif a_file_size is not None and a_file_size == 0:
                raise DataError(f"Input {suffix} file is 0-bytes: {a_file}")
            if a_file_size:
                files.append(a_file)
            else:
                files.append(None)
        return files

    def submit_run(self) -> None:
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {self.model.id}: Submit to Cromwell")
        try:
            infiles = self.get_run_input()
        except DataError as error:
            logger.error(f"Run {self.model.id}: {error}")
            self.update_run_status("submission failed", f"Bad input: {error}")
            return
        try:
            cromwell_run_id = cromwell.submit(*infiles)
        except CromwellError as error:
            logger.error(f"Run {self.model.id} submission failed: {error}")
            self.update_run_status("submission failed", f"{error}")
        else:
            self.model.cromwell_run_id = cromwell_run_id
            self.update_run_status("submitted", f"cromwell_run_id={cromwell_run_id}")

    def check_run_cromwell_status(self) -> None:
        """
        Check Cromwell for the status of the Run.
        """
        logger.debug(f"Run {self.model.id}: Check Cromwell status")
        try:
            cromwell_status = cromwell.get_status(self.model.cromwell_run_id)
        except CromwellError as error:
            logger.error(
                f"Unable to check Cromwell status of Run {self.model.id}: {error}"
            )
            raise

        # check if state has changed; allowed states and transitions for a successful run are:
        # submitted -> queued -> running -> succeeded (with no skipping of states allowed)
        # Additionally, any state can transition directly to cancelled or failed.
        logger.debug(f"Run {self.model.id}: Cromwell status is {cromwell_status}")
        if cromwell_status == "Running":
            # no skips allowed, so there may be more than one transition
            if self.model.status == "submitted":
                self.update_run_status("queued")
            elif self.model.status == "queued":
                # although Cromwell may consider a Run to be "Running", since it does not distinguish between
                # "queued" and "running", we check the task-log to see if any task is "running"; only once any
                # task is running does the Run transition to the "running" state.
                tasks_status = tasks.get_run_status(self.session, self.model.id)
                if tasks_status == "running":
                    self.update_run_status("running")
        elif cromwell_status == "Failed":
            self.update_run_status("failed")
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition
            if self.model.status == "submitted":
                self.update_run_status("queued")
            if self.model.status == "queued":
                self.update_run_status("running")
            self.update_run_status("succeeded")
        elif cromwell_status == "Aborted":
            self.update_run_status("cancelled")

    def _write_outputs_json(self, metadata):
        """Write outputs.json to workflow_root dir"""
        cromwell_workflow_dir = metadata.workflow_root()
        outputs_file = os.path.join(cromwell_workflow_dir, "outputs.json")
        outputs = metadata.outputs(relpath=True)
        with open(outputs_file, "w") as fh:
            fh.write(json.dumps(outputs, sort_keys=True, indent=4))

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

        try:
            self._write_outputs_json(metadata)
        except OSError as error:
            logger.error(f"Run {self.model.id}: Cannot write outputs json: {error}")
            # don't change state; keep trying
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


def send_run_status_logs(session, central_rpc_client) -> None:
    """Send run logs to Central"""

    # get updates from datbase
    try:
        query = (
            session.query(models.Run_Log).filter(models.Run_Log.sent.is_(False)).all()
        )
    except SQLAlchemyError as error:
        logger.exception(f"Unable to select from run_logs: {error}")
        return
    num_logs = len(query)
    if not num_logs:
        return
    logger.debug(f"Sending {num_logs} run logs")

    for log in query:
        data = {
            "site_id": config.conf.get("SITE", "id"),
            "run_id": log.run_id,
            "status_from": log.status_from,
            "status_to": log.status_to,
            "timestamp": log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            "reason": log.reason,
        }
        # add special fields
        if log.status_to == "submitted":
            run = session.query(models.Run).get(log.run_id)
            data["cromwell_run_id"] = run.cromwell_run_id
        elif log.status_to == "downloading":
            run = session.query(models.Run).get(log.run_id)
            data["download_task_id"] = run.download_task_id
        try:
            response = central_rpc_client.request("update_run_logs", data)
        except Exception as error:
            logger.exception(f"RPC update_run_logs error: {error}")
            continue
        if "error" in response:
            logger.info(f"RPC update_run_status failed: {response['error']['message']}")
            continue
        log.sent = True
        try:
            session.commit()
        except Exception as error:
            session.rollback()
            logger.exception(f"Error updating run_logs as sent: {error}")


def get_run_status_logs(session, run_id) -> None:
    """Get run logs for a given run."""

    # get updates from datbase
    try:
        query = session.query(models.Run_Log) \
            .filter(models.Run_Log.run_id == run_id) \
            .all()
    except SQLAlchemyError as error:
        logger.exception(f"Unable to select from run_logs: {error}")
        return

    datas = []
    for log in query:
        data = {
            "site_id": config.conf.get("SITE", "id"),
            "run_id": log.run_id,
            "status_from": log.status_from,
            "status_to": log.status_to,
            "timestamp": log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            "reason": log.reason,
        }
        datas.append(data)
    return datas
