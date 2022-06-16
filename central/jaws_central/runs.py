import logging
import time
from datetime import datetime
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
import json
import smtplib
import ssl
import requests
from jaws_central import models
from jaws_central import config
from jaws_central import jaws_constants
from jaws_central.transfers import Transfer, TransferNotFoundError

logger = logging.getLogger(__package__)


class RunError(Exception):
    # base class for all errors in this package
    pass


class RunDbError(RunError):
    pass


class RunNotFound(RunError):
    pass


class RunSiteError(RunError):
    pass


class DataError(RunError):
    pass


class RunRpcError(RunError):
    pass


class Run:
    """Class representing a single Run"""

    def __init__(self, session, data, rpc_index=None):
        self.session = session
        self.data = data
        self.operations = {
            "created": self.submit_upload,
            "upload queued": self.check_if_upload_complete,
            "uploading": self.check_if_upload_complete,
            "upload complete": self.submit_run,
            "finished": self.submit_download,
            "cancelled": self.submit_download,
            "download queued": self.check_if_download_complete,
            "downloading": self.check_if_download_complete,
            "download complete": self.send_email,
            "email sent": self.post_to_webhook,
        }
        self.rpc_index = rpc_index
        self.upload = None
        self.download = None

    @classmethod
    def from_params(cls, session, **kwargs):
        """Insert run record into rdb"""
        manifest_json = "[]"
        if "manifest" in kwargs:
            # if a list is provided then need to convert to JSON text
            assert type(kwargs["manifest"]) == list
            manifest_json = json.dumps(kwargs["manifest"])
        elif "manifest_json" in kwargs:
            assert type(kwargs["manifest_json"]) == str
            manifest_json = kwargs["manifest_json"]
        try:
            data = models.Run(
                user_id=kwargs["user_id"],
                submission_id=kwargs["submission_id"],
                max_ram_gb=int(kwargs["max_ram_gb"]),
                caching=kwargs["caching"],
                input_site_id=kwargs["input_site_id"],
                compute_site_id=kwargs["compute_site_id"],
                status="uploading",
                wdl_file=kwargs["wdl_file"],
                json_file=kwargs["json_file"],
                tag=kwargs["tag"],
                manifest_json=manifest_json,
                webhook=kwargs["webhook"],
            )
        except SQLAlchemyError as error:
            raise (f"Error creating model for new Run {kwargs['run_id']}: {error}")
        try:
            session.add(data)
            session.commit()
            # data.id will be defined after the commit
        except SQLAlchemyError as error:
            session.rollback()
            raise (error)
        else:
            return cls(session, data)

    @classmethod
    def from_id(cls, session, run_id):
        """Select run record from rdb"""
        try:
            data = session.query(models.Run).get(run_id)
        except IntegrityError as error:
            logger.error(f"Run {run_id} not found: {error}")
            raise RunNotFound(f"Run {run_id} not found")
        except SQLAlchemyError as error:
            err_msg = f"Unable to select run, {run_id}: {error}"
            logger.error(err_msg)
            raise RunDbError(err_msg)
        else:
            return cls(session, data)

    @property
    def status(self) -> str:
        """Return the current state of the run."""
        return self.data.status

    def info(self, verbose: bool = False):
        """
        Return dictionary of run info.
        Run data cannot be changed by altering the returned dict.
        :param verbose: True if more fields desired else fewer.
        :type verbose: bool
        :return: selected fields
        :rtype: dict
        """
        info = {
            "id": self.data.id,
            "compute_site_id": self.data.compute_site_id,
            "result": self.data.result,
            "status": self.data.status,
            "updated": self.data.updated.strftime("%Y-%m-%d %H:%M:%S"),
        }
        if verbose:
            more_info = {
                "submitted": self.data.submitted.strftime("%Y-%m-%d %H:%M:%S"),
                "status_detail": jaws_constants.run_status_msg.get(
                    self.data.status, ""
                ),
                "input_site_id": self.data.input_site_id,
                "cromwell_run_id": self.data.cromwell_run_id,
                "upload_id": self.data.upload_id,
                "submission_id": self.data.submission_id,
                "download_id": self.data.download_id,
                "user_id": self.data.user_id,
                "tag": self.data.tag,
                "wdl_file": self.data.wdl_file,
                "json_file": self.data.json_file,
            }
            info.update(more_info)
        return info

    def check_status(self) -> None:
        """Check the run's status, promote to next state if ready"""
        status = self.data.status
        logger.debug(f"Run {self.data.id} is {status}")
        if status in self.operations:
            return self.operations[status]()

    def cancel(self) -> None:
        """Cancel a run, aborting Cromwell or file transfer as appropriate"""
        if self.data.status == "upload queued":
            upload = self.get_upload()
            upload.cancel()
        elif self.data.status == "download queued":
            download = self.get_download()
            download.cancel()
        elif self.data.status in [
            "submitted",
            "queued",
            "running",
        ]:
            rpc_client = self.rpc_index.get_client(self.data.compute_site_id)
            params = {"user_id": self.data.user_id, "run_id": self.data.id}
            rpc_client.request("cancel", params)
        self.update_status("cancelled")

    def _get_transfer(self, transfer_id: int):
        """
        Init transfer object (either an upload or a download).
        """
        try:
            transfer = Transfer.from_id(self.session, transfer_id)
        except TransferNotFoundError:
            raise
        except Exception:
            raise
        else:
            return transfer

    def get_upload(self):
        """Get transfer object for the upload"""
        if self.data.upload_id is None:
            self.upload = None
        elif self.upload is None:
            self.upload = self._get_transfer(self.data.upload_id)
        return self.upload

    def get_download(self):
        """Get transfer object for the download"""
        if self.data.download_id is None:
            self.download = None
        elif self.download is None:
            self.download = self._get_transfer(self.data.download_id)
        return self.download

    def inputs_manifest(self) -> list:
        """
        The manifest stored in the object and db are just relative paths.  We need the
        absolute paths at the source and destination sites to transfer the files.
        :return: manifest table (each row is src, dest absolute paths)
        :rtype: list
        """
        return json.loads(self.data.manifest_json)

    def outputs_manifest(self) -> list:
        """
        Generate table of absolute paths of the files to transfer.
        :return: manifest table (each row is src, dest absolute paths)
        :rtype: list
        """
        # request list of output files from the compute-site's Cromwell
        rpc_client = self.rpc_index.get_client(self.data.compute_site_id)
        params = {
            "run_id": self.data.id,
            "complete": True,
        }
        try:
            response = rpc_client.request("run_manifest", params)
        except Exception as error:
            logger.error(f"RPC request output manifest for Run {self.data.id}: {error}")
            raise RunRpcError(error)
        if "error" in response:
            error = f"RPC request output manifest for Run {self.data.id}: {response['error']['message']}"
            logger.error(error)
            raise RunRpcError(error)
        result = response["result"]
        manifest = result["manifest"]
        workflow_root = result["workflow_root"]
        return workflow_root, manifest

    def submit_upload(self):
        logger.debug(f"Run {self.data.id}: Submit upload")

        # if input and compute site are same, there are no files to transfer, so just
        # promote the state (two updates are required since we don't skip states)
        if self.data.input_site_id == self.data.compute_site_id:
            self.update_status("upload queued")
            time.sleep(1)
            self.update_status("upload complete")
            return

        src_config = config.conf.get_site(self.data.input_site_id)
        dest_config = config.conf.get_site(self.data.compute_site_id)
        params = {
            "src_site_id": self.data.input_site_id,
            "src_base_dir": src_config.get("inputs_dir"),
            "dest_site_id": self.data.compute_site_id,
            "dest_base_dir": dest_config.get("inputs_dir"),
            "manifest": self.inputs_manifest(),
        }
        try:
            transfer = Transfer.from_params(self.session, params)
            transfer.submit_transfer()
        except Exception as error:
            logger.error(f"Failed to create upload: {error}")
            self.update_status("upload failed", f"{error}")
        else:
            logger.debug(f"Run {self.data.id} upload {transfer.data.id} queued")
            self.data.upload_id = transfer.data.id
            self.update_status("upload queued")

    def submit_download(self):
        logger.debug(f"Run {self.data.id}: Submit download")

        if self.data.cromwell_run_id is None:
            self.update_status(
                "download complete", "The run was cancelled before Cromwell; no output was created."
            )
            return

        # if input and compute site are same, there are no files to transfer, so just
        # promote the state (two updates are required since we don't skip states)
        if self.data.input_site_id == self.data.compute_site_id:
            self.update_status("download queued")
            time.sleep(1)  # hack so they don't have the same timestamp
            self.update_status("download complete")
            return

        dest_config = config.conf.get_site(self.data.input_site_id)
        workflow_root, manifest = self.outputs_manifest()
        dest_base_dir = f"{dest_config.get('downloads_dir')}/{self.data.submission_id}"
        params = {
            "src_site_id": self.data.compute_site_id,
            "src_base_dir": workflow_root,
            "dest_site_id": self.data.input_site_id,
            "dest_base_dir": dest_base_dir,
            "manifest": manifest,
        }
        try:
            transfer = Transfer.from_params(self.session, params)
        except Exception as error:
            logger.error(f"Failed to create download: {error}")
            self.update_status("download failed", f"{error}")
            return
        try:
            transfer.submit_transfer()
        except Exception as error:
            logger.error(f"Failed to submit download: {error}")
            self.update_status("download failed", f"{error}")
        else:
            logger.debug(f"Run {self.data.id} download {transfer.data.id} queued")
            self.data.download_id = transfer.data.id
            self.update_status("download queued")

    def check_if_upload_complete(self) -> None:
        """
        If the upload is complete, update the run's status.
        """
        logger.debug(f"Run {self.data.id}: Check upload status")
        try:
            upload = self.get_upload()
        except Exception as error:
            logger.error(f"Unable to get upload for Run {self.data.id}: {error}")
            return
        status = upload.status()
        if status in ["submission failed", "failed"]:
            self.update_status("upload failed")
        elif status == "succeeded":
            self.update_status("upload complete")

    def check_if_download_complete(self) -> None:
        """
        If the download is complete, update the run's status.
        """
        logger.debug(f"Run {self.data.id}: Check download status")
        try:
            download = self.get_download()
        except Exception as error:
            logger.error(f"Unable to get download for Run {self.data.id}: {error}")
            return
        status = download.status()
        if status in ["submission failed", "failed"]:
            self.update_status("download failed")
        elif status == "succeeded":
            time.sleep(
                30
            )  # give file system a chance to update all metadata tables (see: #1236)
            self.update_status("download complete")

    def user_email(self):
        """
        Get user's email from RDb.
        """
        try:
            user = (
                self.session.query(models.User)
                .filter(models.User.id == self.data.user_id)
                .one_or_none()
            )
        except SQLAlchemyError as error:
            logger.warning(f"Failed to select user email: {error}", exc_info=True)
            raise
        else:
            return user.email

    def submit_run(self) -> None:
        """
        After the upload is complete, send the Run to the compute jaws-site via RPC
        """
        site_id = self.data.compute_site_id
        logger.debug(f"Run {self.data.id}: Submit to {site_id}")
        rpc_client = self.rpc_index.get_client(site_id)
        params = {
            "user_id": self.data.user_id,
            "run_id": self.data.id,
            "caching": self.data.caching,
            "submission_id": self.data.submission_id,
            "input_site_id": self.data.input_site_id,
        }
        # submit to compute-site via RPC but do not fail upon error; try again later instead
        try:
            response = rpc_client.request("submit_run", params)
        except Exception as error:
            logger.error(f"RPC submit run to {site_id} failed: {error}")
            return
        if "error" in response:
            error = f"Run {self.data.id} rejected by {site_id}: {response['error']['message']}"
            logger.error(error)
            return
        self.update_status("ready")

    def update_status(self, status_to, reason=None) -> None:
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        status_from = self.data.status
        logger.info(f"Run {self.data.id}: now {status_to}")
        timestamp = datetime.utcnow()
        self._update_status(status_to, timestamp)
        self._insert_run_log(status_from, status_to, timestamp, reason)

    def _update_status(self, new_status, timestamp) -> None:
        """
        Update Run's current status in 'runs' table
        """
        try:
            self.data.status = new_status
            self.data.updated = timestamp
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Unable to update Run {self.data.id}: {error}")
            raise

    def _insert_run_log(self, status_from, status_to, timestamp, reason) -> None:
        """
        Save record of state transition in 'run_logs' table
        """
        try:
            log_entry = models.Run_Log(
                run_id=self.data.id,
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
                f"Failed to insert log for Run {self.data.id} ({status_to}): {error}"
            )
            raise

    def send_email(self):
        receiver_email = self.user_email()

        if receiver_email is None:
            self.update_status("email sent")
            return

        sender_email = config.conf.get("EMAIL", "user")
        smtp_server = config.conf.get("EMAIL", "server")
        port = config.conf.get("EMAIL", "port")
        password = config.conf.get("EMAIL", "password")

        tag_text = f" ({self.data.tag})" if self.data.tag else ""

        message = f"""Subject: JAWS Run {self.data.id} {self.data.result}{tag_text}

        Your run has completed.

        run_id: {self.data.id}
        result: {self.data.result}
        wdl_file: {self.data.wdl_file}
        json_file: {self.data.json_file}
        tag: {self.data.tag}
        """

        context = ssl.create_default_context()
        try:
            with smtplib.SMTP(smtp_server, port) as server:
                server.starttls(context=context)
                server.login(sender_email, password)
                server.sendmail(sender_email, receiver_email, message)
        except Exception as error:
            logger.error(f"Run {self.data.id} failed to send email: {error}")
        else:
            self.update_status("email sent")

    def post_to_webhook(self):
        if self.data.webhook:
            try:
                requests.post(self.data.webhook, data=self.info())
            except Exception as error:
                self.error(f"Run {self.data.id} failed to POST: {error}")
            else:
                self.update_status("done")
        else:
            self.update_status("done")


def check_active_runs(session, rpc_index) -> None:
    """
    Get active runs from db and have each check and update their status.
    """
    active_states = [
        "created",
        "upload queued",
        "uploading",
        "upload complete",
        "finished",
        "cancelled",
        "download queued",
        "downloading",
        "download complete",
        "email sent",
    ]
    try:
        rows = (
            session.query(models.Run).filter(models.Run.status.in_(active_states)).all()
        )
    except SQLAlchemyError as error:
        logger.warning(f"Failed to select active runs from db: {error}", exc_info=True)
    else:
        if len(rows):
            logger.debug(f"Central is responsible for {len(rows)} run(s)")
            for row in rows:
                run = Run(
                    session,
                    row,
                    rpc_index,
                )
                run.check_status()
