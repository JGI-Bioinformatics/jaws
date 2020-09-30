"""
JAWS Analysis Service API
"""

import logging
import os
import requests
from datetime import datetime
import sqlalchemy.exc
from sqlalchemy.exc import SQLAlchemyError
import globus_sdk
from jaws_site import config, db
from jaws_site.cromwell import Cromwell
from jaws_rpc.rpc_client import RpcClient


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    """Generic exception when Cromwell does not return metadata."""

    pass


class RunNotFound(Exception):
    pass


class CromwellError(Exception):
    pass


class Run:
    """
    Representation of a run, which is the execution of a workflow (WDL) on a specific set of inputs (JSON) by Cromwell.
    No other object has direct access to Cromwell.
    A run's persistent data is stored in a relational database.
    If the Run does not already exist, the parameters to initialize it are required.
    """

    client_id = config.conf.get("GLOBUS", "client_id")

    def __init__(self, run_id: int, session=None, params=None):
        self.run_id = run_id
        self.data = None  # row in Runs table
        if session:
            self.session = session
        else:
            self.session = db.Session()
        if params:
            # this is a new Run; insert into database.
            self.__add_run__(params)
        else:
            # retrieve existing Run from db or raise exception
            self.__load_run__()

    def __add_run__(self, params):
        """
        Save new run submission in database.
        New runs are always initialized in the "uploading" state.
        """
        user_id = params["user_id"]
        logger.info(f"User {user_id}: Add Run {self.run_id}")
        try:
            run = db.Run(
                id=self.run_id,
                user_id=params["user_id"],
                submission_id=params["submission_id"],
                upload_task_id=params["upload_task_id"],
                output_endpoint=params["output_endpoint"],
                output_dir=params["output_dir"],
                status="uploading",
            )
        except Exception as error:
            logger.exception(f"Invalid new run parameters, {params}: {error}")
            raise ValueError(f"Invalid Run parameters: {error}")
        try:
            self.session.add(run)
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Failed to insert new Run record: {error}")
            raise DatabaseError(f"Error saving new run: {error}")
        self.data = run

    def __load_run__(self):
        """
        Retrieve the existing Run's record from the database.
        """
        try:
            run = self.session.query(db.Run).get(self.run_id)
        except sqlalchemy.exc.IntegrityError as error:
            logger.exception(f"Run not found: {self.run_id}: {error}")
            raise RunNotFound(f"{error}")
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting on runs table: {error}")
            raise DatabaseError(f"{error}")
        if not run:
            logger.debug(f"Run {self.run_id} not found")
            raise RunNotFound(f"No such record")
        self.data = run

    def get_status(self):
        """
        Return the current status of the Run.
        """
        # TODO EXTRACT FROM RUN LOG INSTEAD AND SIMPLIFY UPDATE_RUN_LOGS METHOD
        return self.data.status

    def get_log(self):
        """
        Get complete Run Log.
        """
        try:
            logs = self.session.query(db.Run_Log).filter_by(run_id=self.run_id).all()
        except sqlalchemy.exc.IntegrityError as error:
            logger.exception(f"Run Logs not found for Run {self.run_id}: {error}")
            raise RunNotFound(f"{error}")
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting on run_logs table: {error}")
            raise DatabaseError(f"{error}")
        if not logs:
            logger.debug(f"Run {self.run_id} not found in logs")
            raise RunNotFound(f"No such record")
        return logs

    def get_metadata(self):
        """Retrieve the metadata of a run.

        :param run_id: Run ID
        :type params: dict
        :return: The Cromwell metadata for the specified run (may be None)
        :rtype: dict
        """
        logger.info(f"Get metadata for Run {self.run_id}")
        if self.data.cromwell_run_id is None:
            return None

        try:
            result = cromwell.get_all_metadata(self.data.cromwell_run_id)
        except requests.exceptions.HTTPError as error:
            logger.exception(f"Get metadata for {self.run_id} failed: {error}")
            raise CromwellError(f"{error}")
        return result

    def cancel(self):
        """Cancel a run.

        :param run_id: JAWS Run ID
        :type run_id: int
        """
        logger.info(f"Cancel run {self.run_id}")

        status = self.data.status
        self.data.status = "cancelled"
        try:
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Error updating Run {self.run_id}: {error}")
            raise DatabaseError(f"{error}")

        # tell Cromwell to cancel the run if it has been submitted to Cromwell already
        if self.data.cromwell_run_id and status in ["submitted", "queued", "running"]:
            logger.debug(f"Run {self.run_id} is {status}: Instructing Cromwell to cancel")
            try:
                cromwell.abort(self.data.cromwell_run_id)
            except requests.exceptions.HTTPError as error:
                logger.exception(
                    f"Error aborting cromwell {self.data.cromwell_run_id}: {error}"
                )
                raise CromwellError(
                    f"The run was cancelled but Cromwell returned an error: {error}"
                )

    def get_errors(self):
        """Retrieve error messages and stderr for failed Tasks.
        :return: error messages and stderr for failed Tasks
        :rtype: dict
        """
        logger.info(f"Get errors for Run {self.run_id}")
        if self.data.cromwell_run_id is None:
            return None
        try:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            result = metadata.errors()
        except requests.exceptions.HTTPError as error:
            logger.exception(f"Get errors for {self.run_id} failed: {error}")
            raise CromwellError(f"{error}")
        return result

    def get_user_service_client(self):
        if not self.user_service:
            self.user_service = RpcClient(config.conf.user_service_params())
        return self.user_service

    def get_transfer_client(self):
        """
        Get a Globus transfer client, authorized to transfer files on the user's behalf.
        """
        if not self.transfer_client:
            # get Globus token from JAWS User service
            user_id = self.get_user_id()
            user_service = self.get_user_service()
            response = user_service.call("get_transfer_token", {"user_id": user_id})
            token = response["result"]

            # get Globus transfer client
            client = globus_sdk.NativeAppAuthClient(self.client_id)
            authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
            self.transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
        return self.transfer_client

    def _get_globus_transfer_status(self, task_id):
        """
        Query Globus transfer service for transfer task status.
        """
        transfer_client = self.get_transfer_client()
        try:
            task = transfer_client.get_task(task_id)
            globus_status = task["status"]
        except Exception:
            logger.exception("Failed to check Globus upload status", exc_info=True)
            raise
        return globus_status

    def check_if_upload_complete(self):
        """
        Check if Globus upload task is complete.
        """
        logger.debug(f"Run {self.run_id}: Check upload status")
        upload_task_id = self.get_upload_task_id()
        try:
            globus_status = self._get_globus_transfer_status(upload_task_id)
        except Exception as error:
            logger.exception(
                f"Failed to check upload {self.data.upload_task_id}: {error}"
            )
            return
        if globus_status == "FAILED":
            self.update_run_status("upload failed")
        elif globus_status == "INACTIVE":
            self.update_run_status(
                "upload inactive", "Your endpoint authorization has expired"
            )
        elif globus_status == "SUCCEEDED":
            self.update_run_status("upload complete")

    def submit_run(self, run):
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {self.run_id}: Submit to Cromwell")

        # Validate input
        file_path = os.path.join(
            self.uploads_dir, self.data.user_id, self.data.submission_id
        )
        wdl_file = file_path + ".wdl"
        json_file = file_path + ".json"
        zip_file = file_path + ".zip"  # might not exist
        if not os.path.exists(json_file):
            logger.warning(f"Missing inputs JSON for run {self.run_id}: {json_file}")
            self.update_run_status("missing input", "Missing inputs JSON file")
            return
        if not os.path.exists(wdl_file):
            logger.warning(f"Missing WDL for run {self.run_id}: {wdl_file}")
            self.update_run_status("missing input", "Missing WDL file")
            return
        if not os.path.exists(zip_file):
            zip_file = None

        # submit to Cromwell
        cromwell_run_id = cromwell.submit(wdl_file, json_file, zip_file)
        if cromwell_run_id:
            self.data.cromwell_run_id = cromwell_run_id
            self.update_run_status(
                "submitted", f"cromwell_run_id={self.data.cromwell_run_id}"
            )
        else:
            self.update_run_status("submission failed")

    def check_run_cromwell_status(self, run):
        """
        Check Cromwell for the status of one Run.
        """
        logger.debug(f"Run {self.run_id}: Check Cromwell status")
        cromwell_status = None
        if self.data.cromwell_workflow_dir:
            # all we need is status
            try:
                cromwell_status = cromwell.get_status(self.data.cromwell_run_id)
            except Exception as error:
                logger.error(
                    f"Unable to check Cromwell status of Run {self.run_id}: {error}"
                )
                return  # try again next time
        else:
            # get complete metadata
            try:
                metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            except Exception:
                return  # try again next time
            if metadata is None:
                return
            self.data.cromwell_workflow_dir = metadata.get("workflowRoot")
            self.session.commit()
            cromwell_status = metadata.get("status")

        # check if state has changed; allowed states and transitions for a successful run are:
        # submitted -> queued -> running -> succeeded (with no skipping of states allowed)
        # Additionally, any state can transition directly to cancelled or failed.
        logger.debug(f"Run {self.run_id}: Cromwell status is {cromwell_status}")
        if cromwell_status == "Running":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.update_run_status("queued")
        elif cromwell_status == "Failed":
            self.update_run_status("failed")
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.update_run_status("queued")
            if self.data.status == "queued":
                self.update_run_status("running")
            if self.data.status == "running":
                self.update_run_status("succeeded")
        elif cromwell_status == "Aborted":
            self.update_run_status("cancelled")

    def __transfer_folder(self, label, transfer_rt, src_dir, dest_endpoint, dest_dir):
        """
        Recursively transfer folder via Globus
        :param label: Label to attach to transfer (e.g. "Run 99")
        :type label: str
        :param transfer_rt: User's Globus transfer refresh token
        :type transfer_rt: str
        :param src_dir: Folder to transfer
        :type src_dir: str
        :param dest_endpoint: Globus endpoint for destination
        :type dest_endpoint: str
        :param dest_dir: Destination path
        :type dest_dir: str
        :return: Globus transfer task id
        :rtype: str
        """
        logger.debug(f"Globus xfer {label}")
        if not src_dir.startswith(self.globus_root_dir):
            logger.error(f"Dir is not accessible via Globus: {src_dir}")
            return None
        rel_src_dir = os.path.relpath(src_dir, self.globus_default_dir)
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to get Globus transfer client to xfer {label}", exc_info=True
            )
            return None
        try:
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.globus_endpoint,
                dest_endpoint,
                label=label,
                sync_level="exists",
                verify_checksum=False,
                preserve_timestamp=True,
                notify_on_succeeded=False,
                notify_on_failed=False,
                notify_on_inactive=False,
                skip_activation_check=True,
            )
            if self.globus_root_dir == "/":
                tdata.add_item(src_dir, dest_dir, recursive=True)
            else:
                tdata.add_item(rel_src_dir, dest_dir, recursive=True)
        except Exception:
            logger.warning(
                f"Failed to prepare download manifest for {label}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to download results with Globus for {label}", exc_info=True,
            )
            return None
        return transfer_result["task_id"]

    def transfer_results(self, run, metadata=None):
        """
        Send run output via Globus
        """
        logger.debug(f"Run {self.run_id}: Download output")
        if self.data.cromwell_workflow_dir is None:
            if not metadata:
                try:
                    metadata = cromwell.get_metadata(run.cromwell_run_id)
                except Exception:
                    return
            if metadata is None:
                return
            outdir = metadata.get("workflowRoot")
            if outdir:
                self.data.cromwell_workflow_dir = outdir
                self.session.commit()
            else:
                # This run failed before a folder was created; nothing to xfer
                self.update_run_status("download complete", "No run folder was created")
                return
        transfer_rt = self.get_transfer_refresh_token
        transfer_task_id = self.__transfer_folder(
            f"Run {self.run_id}",
            transfer_rt,
            self.data.cromwell_workflow_dir,
            self.data.output_endpoint,
            self.data.output_dir,
        )
        if transfer_task_id:
            self.data.download_task_id = transfer_task_id
            self.update_run_status(
                "downloading", f"download_task_id={self.data.download_task_id}"
            )
            self.session.commit()
        # else ignore error; was logged; will try again later

    def check_if_download_complete(self, run):
        """
        If download is complete, change state.
        """
        logger.debug(f"Run {self.run_id}: Check download status")
        try:
            globus_status = self._get_globus_transfer_status(
                run, self.data.download_task_id
            )
        except Exception as error:
            logger.exception(
                f"Failed to check download {self.data.download_task_id}: {error}"
            )
            return
        if globus_status == "SUCCEEDED":
            self.update_run_status("download complete")
        elif globus_status == "FAILED":
            self.update_run_status("download failed")

    def update_run_status(self, new_status, reason=None):
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        logger.info(f"Run {self.run_id}: now {new_status}")
        status_from = self.data.status
        timestamp = datetime.utcnow()

        # update current status in "runs" table
        try:
            self.data.status = new_status
            self.data.updated = timestamp
            self.session.commit()
        except Exception as error:
            logger.exception(f"Unable to update Run {self.run_id}: {error}")

        # populate "reason" field
        if new_status == "submitted":
            reason = f"cromwell_run_id={self.data.cromwell_run_id}"
        elif new_status == "downloading":
            reason = f"download_task_id={self.data.download_task_id}"

        # save state transition in "run_logs" table
        try:
            log_entry = db.Run_Log(
                run_id=self.run_id,
                status_from=status_from,
                status_to=new_status,
                timestamp=timestamp,
                reason=reason,
            )
            self.session.add(log_entry)
            self.session.commit()
        except Exception as error:
            logger.exception(
                f"Failed to insert log for Run {self.run_id} : {new_status}, {reason}: {error}"
            )
        # notifying Central of state change is handled by send_run_status_logs
