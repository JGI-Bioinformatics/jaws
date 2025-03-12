# Class representing Runs
#
# Description of Possible States:
# - ready* : Initial state.  Infiles have been already been uploaded to this Site.
# - submission failed : The Run could not be submitted to Cromwell.
# - submitted : The Run was successfully submitted to Cromwell and a cromwell_run_id was returned.
# - queued : No tasks for the Run have started executing yet.
# - running : At least one task of the Run has started running.
# - failed : The Run was failed by Cromwell.
# - succeeded : The Run has completed successfully by Cromwell.
# - complete : Supplementary files have been added to the Run's workflow root dir
# - finished** : The Run metrics report has been published to ElasticSearch (if wasn't cancelled).
# - cancel : Mark a Run to be cancelled (via RPC from Central); will be cancelled by run-daemon later.
# - cancelled** : The Run has been successfully cancelled.
#
# * this is the only initial state
# ** these are the two terminal states


import io
import json
import logging
import os
import shutil
from datetime import datetime
from random import shuffle
import boto3
import botocore
from tenacity import retry, stop_after_attempt, wait_fixed, before_log, after_log

from jaws_site import config, models
from jaws_site.cromwell import (
    Cromwell,
    CromwellError,
    CromwellRunError,
    CromwellRunNotFoundError,
    CromwellServiceError,
    CromwellGetMetadataError,
)
from jaws_site.tasks import TaskLog
from jaws_site.utils import write_json_file
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.session import sessionmaker

logger = logging.getLogger(__package__)


if config.conf is not None:
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
else:
    cromwell = Cromwell("localhost")

MAX_ERROR_STRING_LEN = 1024
JAWS_GET_METADATA_MAX_RETRIALS = int(
    config.conf.get("SITE", "jaws_get_metadata_max_retrials", 3)
)
JAWS_GET_METADATA_WAIT_SEC = int(
    config.conf.get("SITE", "jaws_get_metadata_wait_sec", 180)
)


def set_atime_now(path: str) -> None:
    """
    Set the file's access time to now.
    """
    if not os.path.isfile(path):
        raise IOError(f"File not found: {path}")
    stat = os.stat(path)
    os.utime(path, times=(datetime.now().timestamp(), stat.st_mtime))


class RunDbError(Exception):
    pass


class RunNotFoundError(Exception):
    pass


class RunFileNotFoundError(Exception):
    pass


class RunInputError(Exception):
    pass


class Run:
    """Class representing a single Run"""

    def __init__(self, session, data, **kwargs):
        self.session = session
        self.data = data
        self._metadata = None
        self.operations = {
            "ready": self.submit_run,
            "submitted": self.check_cromwell_run_status,
            "queued": self.check_cromwell_run_status,
            "running": self.check_cromwell_run_status,
            "succeeded": self.write_supplement,
            "failed": self.write_supplement,
            "cancelled": self.write_supplement,
            "complete": self.publish_report,
            "cancel": self.cancel,
        }
        self.central_rpc_client = kwargs.get("central_rpc_client")
        self.reports_rpc_client = kwargs.get("reports_rpc_client")

        try:
            self.config = {
                "site_id": kwargs.get("site_id", config.conf.get("SITE", "id")),
                "inputs_dir": kwargs.get("inputs_dir", config.conf.get("SITE", "inputs_dir")),
                "default_container": kwargs.get("default_container", config.conf.get(
                    "SITE", "default_container", "ubuntu:latest"
                )),
                "max_user_active_runs": kwargs.get("max_user_active_runs", int(
                    config.conf.get("SITE", "max_user_active_runs", 0)
                )),
                "aws_access_key_id": kwargs.get("aws_access_key_id", config.conf.get("AWS", "aws_access_key_id")),
                "aws_region_name": kwargs.get("aws_region_name", config.conf.get("AWS", "aws_region_name")),
                "aws_secret_access_key": kwargs.get("aws_secret_access_key", config.conf.get(
                    "AWS", "aws_secret_access_key"
                )),
                "cromwell_url": kwargs.get("cromwell_url", config.conf.get("CROMWELL", "url")),
                "cromwell_executions_dir": kwargs.get("cromwell_executions_dir", config.conf.get(
                    "CROMWELL", "executions_dir"
                )),
            }
        except Exception as error:
            logger.error(f"Error loading config: {error}")

    @classmethod
    def from_params(
        cls, session, params, central_rpc_client=None, reports_rpc_client=None
    ):
        """Insert new Run into RDb.  Site only receives Runs in the "upload complete" state."""
        try:
            if not isinstance(params["caching"], bool):
                raise SQLAlchemyError
            data = models.Run(
                id=int(params["run_id"]),
                user_id=params["user_id"],
                caching=params["caching"],
                submission_id=params["submission_id"],
                input_site_id=params["input_site_id"],
                status="ready",
                wdl_basename=params["wdl_basename"],
                json_basename=params["json_basename"],
                tag=params.get("tag", None),
            )
        except SQLAlchemyError as error:
            raise RunDbError(
                f"Error creating model for new Run {params['run_id']}: {error}"
            )
        try:
            session.add(data)
            session.commit()
        except SQLAlchemyError as error:
            session.rollback()
            raise RunDbError(error)
        else:
            return cls(
                session,
                data,
                central_rpc_client=central_rpc_client,
                reports_rpc_client=reports_rpc_client,
            )

    @classmethod
    def from_id(cls, session, run_id, central_rpc_client=None, reports_rpc_client=None):
        """Select Run record from RDb given primary key"""
        try:
            data = session.query(models.Run).get(run_id)
        except IntegrityError as error:
            logger.warning(f"Run {run_id} not found: {error}")
            raise RunNotFoundError(f"Run {run_id} not found")
        except SQLAlchemyError as error:
            err_msg = f"Unable to select run, {run_id}: {error}"
            logger.error(err_msg)
            raise RunDbError(err_msg)
        else:
            return cls(
                session,
                data,
                central_rpc_client=central_rpc_client,
                reports_rpc_client=reports_rpc_client,
            )

    @classmethod
    def from_cromwell_run_id(
        cls, session, cromwell_run_id, central_rpc_client=None, reports_rpc_client=None
    ):
        """Select Run record from RDb given Cromwell Run ID"""
        try:
            data = (
                session.query(models.Run)
                .filter(models.Run.cromwell_run_id == cromwell_run_id)
                .one()
            )
        except NoResultFound as error:
            logger.warning(f"Cromwell {cromwell_run_id} not found: {error}")
            raise RunNotFoundError(f"Cromwell {cromwell_run_id} not found")
        except SQLAlchemyError as error:
            err_msg = f"Unable to select cromwell, {cromwell_run_id}: {error}"
            logger.error(err_msg)
            raise RunDbError(err_msg)
        else:
            return cls(
                session,
                data,
                central_rpc_client=central_rpc_client,
                reports_rpc_client=reports_rpc_client,
            )

    def status(self) -> str:
        """Return the current state of the run."""
        status = self.data.status
        return status

    def summary(self) -> dict:
        """Produce summary of Run info"""
        summary = {
            "run_id": self.data.id,
            "user_id": self.data.user_id,
            "cromwell_run_id": self.data.cromwell_run_id,
            "json_basename": self.data.json_basename,
            "tag": self.data.tag,
            "workflow_root": self.data.workflow_root,
            "workflow_name": self.data.workflow_name,
            "submitted": self.data.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": self.data.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "status": self.data.status,
            "result": self.data.result,
            "compute_site_id": self.config["site_id"],
            "cpu_hours": self.data.cpu_hours,
        }
        return summary

    def did_run_start(self):
        """
        Check if any task has started running by checking the task log.
        """
        task_log = self.task_log()
        if task_log:
            return task_log.did_run_start()
        else:
            return False

    def task_log(self):
        if self.data.cromwell_run_id is None:
            return None
        else:
            return TaskLog(self.session, self.data.cromwell_run_id)

    def check_status(self) -> None:
        """Check the run's status, promote to next state if ready"""
        status = self.data.status
        if status in self.operations:
            return self.operations[status]()

    def mark_to_cancel(self) -> None:
        """
        Tag a run to be cancelled.  Raise if not successful.
        """
        if self.data.status in ["cancel", "cancelled"]:
            return
        elif self.data.status == "ready":
            # this run wasn't submitted to Cromwell yet, so can be cancelled immediately
            self.update_run_status("cancel")
            self.update_run_status("cancelled")
        elif self.data.cromwell_run_id and self.data.status in [
            "submitted",
            "queued",
            "running",
        ]:
            # an active Cromwell run requires communication with Cromwell, so just mark to cancel for now.
            try:
                self.update_run_status("cancel")
            except Exception as error:
                logger.error(
                    f"Failed to mark Run {self.data.id} to be cancelled: {error}"
                )
                raise RunDbError(
                    f"Change Run {self.data.id} status to 'cancel' failed: {error}"
                )
        else:
            # if the run has already completed, it cannot be cancelled.
            raise RunInputError(
                f"Run {self.data.id}: cannot abort run in {self.data.status} state"
            )

    def cancel(self) -> None:
        """
        Cancel a run, aborting Cromwell if appropriate.

        Note: `cancel` = to-be-cancelled status
        """
        try:
            logger.debug(f"Run {self.data.id}: Send Cromwell abort")
            result = cromwell.abort(self.data.cromwell_run_id)
        except CromwellRunNotFoundError as error:
            logger.error(
                f"Cromwell couldn't cancel unknown Run {self.data.id}: {error}"
            )
            # do not return; future attempts will similarly fail; proceed to mark as cancelled
        except CromwellServiceError as error:
            logger.warning(
                f"Could not cancel Run {self.data.id} as Cromwell unavailable: {error}"
            )
            return  # try again later
        except CromwellError as error:
            logger.error(
                f"Unknown Cromwell error cancelling Run {self.data.id}: {error}"
            )
            # unknown error; don't bother trying again later, just proceed to mark as cancelled
        else:
            logger.debug(f"Run {self.data.id}: Cromwell abort successful: {result}")

        try:
            self.update_run_status("cancelled")
        except Exception as error:
            logger.error(f"Failed to cancel Run {self.data.id}: {error}")
            raise RunDbError(
                f"Change Run {self.data.id} status to cancelled failed: {error}"
            )

    @staticmethod
    def s3_parse_path(full_path):
        full_path = full_path.replace("s3://", "", 1)  # discard S3 identifier
        folders = full_path.split("/")
        s3_bucket = folders.pop(0)
        path = "/".join(folders)
        return s3_bucket, path

    def _read_file_s3(self, path, binary=False):
        """
        Read the contents of a file from S3 into RAM and return a file handle object.
        """
        s3_bucket, src_path = self.s3_parse_path(path)
        logger.debug(f"Read from s3://{s3_bucket} obj {src_path}")

        aws_session = boto3.Session(
            aws_access_key_id=self.config["aws_access_key_id"],
            aws_secret_access_key=self.config["aws_secret_access_key"],
            region_name=self.config["aws_region_name"],
        )
        s3_resource = aws_session.resource("s3")
        bucket_obj = s3_resource.Bucket(s3_bucket)
        fh = io.BytesIO()
        try:
            bucket_obj.download_fileobj(src_path, fh)
        except botocore.exceptions.ClientError as error:
            raise RunFileNotFoundError(f"File obj not found, {src_path}: {error}")
        except Exception as error:
            raise RunFileNotFoundError(error)
        fh.seek(0)

        if not binary:
            data = fh.read()
            fh.close()
            fh = io.StringIO(data.decode("utf-8"))
        return fh

    def _read_file_nfs(self, path: str, binary=False):
        """
        Read file from NFS into RAM and return a file handle object.
        """
        if not os.path.isfile(path):
            raise RunFileNotFoundError(f"File not found: {path}")
        data = None
        mode = "rb" if binary else "r"
        try:
            with open(path, mode) as fh:
                data = fh.read()
        except RunFileNotFoundError:
            raise
        if len(data) == 0:
            raise RunFileNotFoundError("File is 0 bytes")
        fh = io.BytesIO(data) if binary else io.StringIO(data)
        fh.seek(0)
        return fh

    def _read_file(self, path: str, binary=False):
        """
        Read file from NFS or S3 and return contents.
        """
        if self.config["inputs_dir"].startswith("s3://"):
            return self._read_file_s3(path, binary)
        else:
            try:
                return self._read_file_nfs(path, binary)
            except PermissionError as error:
                logger.error(f"Run {self.data.id}: Unable to read {path}: {error}")
                raise

    def _read_json_file(self, path) -> dict:
        """
        Read JSON file and return decoded contents as variable.
        """
        return json.load(self._read_file(path))

    # TODO SWITCH FROM KWARGS TO CONFIG
    def _write_file_s3(self, path: str, content: str, **kwargs):
        s3_bucket, src_path = s3_parse_uri(path)
        aws_session = boto3.Session(
            aws_access_key_id=kwargs["aws_access_key_id"],
            aws_secret_access_key=kwargs["aws_secret_access_key"],
            region_name=kwargs["aws_region_name"],
        )
        s3_resource = aws_session.resource("s3")
        bucket_obj = s3_resource.Bucket(s3_bucket)
        bucket_obj.put(Body=content)

    @staticmethod
    def _write_file_nfs(path: str, content: str):
        try:
            with open(path, "w") as fh:
                fh.write(content)
        except PermissionError as error:
            logger.error(f"Could not write to {path}: {error}")
            raise

    def _write_file(self, path: str, contents: str, **kwargs):
        """
        Write contents to NFS or S3 file.
        :param path: Path to file (may be s3 item)
        :ptype path: str
        :param contents: Contents to write
        :ptype contents: str
        """
        if path.startswith("s3://"):
            return self._write_file_s3(path, contents, **kwargs)
        else:
            return self._write_file_nfs(path, contents)

    def read_inputs(self):
        """
        Read inputs json from file or S3 bucket and return contents.
        :return: The Run's input parameters
        :rtype: dict
        """
        fh = self._read_file(
            os.path.join(self.config["inputs_dir"], f"{self.data.submission_id}.json")
        )
        try:
            inputs = json.load(fh)
        except Exception as error:
            raise ValueError(f"Invalid inputs JSON: {error}")
        return inputs

    def rel_to_abs(self, data: any, root: str) -> any:
        """
        Recursively traverse data structure and replace FILE variables' relative paths
        (indicated by starting with "./") to absolute paths, using this site's inputs dir.
        Also change the access times to ensure the files are not purged.
        Raise on missing file.
        :param prefix: The path to the input data folder
        :ptype prefix: str
        :return: If item is a relpath, return abspath, else return unmodified item.
        :rtype: any
        """
        if type(data) is str:
            if data.startswith("./"):
                abspath = os.path.normpath(os.path.join(root, data))
                if not os.path.isfile(abspath):
                    raise RunFileNotFoundError(f"{data} looks like a file path but was not defined with \"File\" in the WDL. Either change to \"File\" or remove the leading \"./\".")
                try:
                    set_atime_now(abspath)
                except PermissionError as error:
                    logger.error(
                        f"Run {self.data.id}: Unable to change atime for {abspath}: {error}"
                    )
                return abspath
            else:
                return data
        elif type(data) is list:
            new_data = []
            for item in data:
                new_data.append(self.rel_to_abs(item, root))
            return new_data
        elif type(data) is dict:
            new_data = {}
            for key, value in data.items():
                new_key = self.rel_to_abs(key, root)
                new_value = self.rel_to_abs(value, root)
                new_data[new_key] = new_value
            return new_data
        else:
            return data

    def inputs(self):
        """
        Get the Run's relpath inputs and convert to valid abspaths for this Site.
        :return: input parameters
        :rtype: dict
        """
        inputs = self.read_inputs()
        return self.rel_to_abs(inputs, self.config["inputs_dir"])

    def inputs_fh(self):
        """
        Get valid inputs json for this Site and create a file handle, which is
        required to POST it to the Cromwell server (instead of writing a tmpfile).
        """
        inputs = self.inputs()
        json_str = json.dumps(inputs)
        fh = io.StringIO(json_str)
        fh.seek(0)
        return fh

    def get_run_inputs(self):
        """
        Get file handles.
        """
        file_handles = {}
        file_handles["inputs"] = self.inputs_fh()
        path = os.path.join(self.config["inputs_dir"], f"{self.data.submission_id}.wdl")
        file_handles["wdl"] = self._read_file(path)
        try:
            path = os.path.join(
                self.config["inputs_dir"], f"{self.data.submission_id}.zip"
            )
            sub = self._read_file(path, True)
        except Exception:
            logger.debug(f"Run {self.data.id}: Zipped subworkflows not found: {path}")
            pass  # subworkflows are optional
        else:
            file_handles["subworkflows"] = sub
        return file_handles

    def cromwell_options(self):
        default_container = self.config["default_container"]
        options = {
            "caching": self.data.caching,
            "default_container": default_container,
            "hogGroup": self.data.user_id,
        }
        return options

    def submit_run(self) -> None:
        """
        Submit a run to Cromwell.
        """

        logger.debug(f"Run {self.data.id}: Submit to Cromwell")

        # If this user has too many concurrent runs, don't submit any more at this time.
        # This is to enforce fair-sharing; the limit is a configuration parameter.
        if max_active_runs_exceeded(
            self.session, self.data.user_id, self.config["max_user_active_runs"]
        ):
            logger.debug(
                f"Run {self.data.id}: User reached concurrent runs limit; skipping."
            )
            return

        try:
            file_handles = self.get_run_inputs()
        except ValueError as error:
            logger.error(
                f"Run {self.data.id}: Submission failed; invalid inputs JSON: {error}"
            )
            self.update_run_status("submission failed", f"Invalid inputs JSON: {error}")
            return
        except RunFileNotFoundError as error:
            logger.error(f"Run {self.data.id}: Submission failed; {error}")
            self.update_run_status("submission failed", str(error))
            return
        try:
            options = self.cromwell_options()
        except Exception as error:
            logger.error(f"Run {self.data.id}: Submission failed; {error}")
            self.update_run_status("submission failed", f"Options error: {error}")
            return
        try:
            cromwell_run_id = cromwell.submit(file_handles, options)
        except CromwellRunError as error:
            logger.error(
                f"Run {self.data.id}: Submission failed (value error): {error}"
            )
            self.update_run_status("submission failed", f"{error}")
        except CromwellError as error:
            logger.error(
                f"Run {self.data.id}: Submission failed (service unavailable): {error}"
            )
            # do not change status to failed; try again later
            return
        else:
            self.data.cromwell_run_id = cromwell_run_id
            self.update_run_status("submitted")

    def resubmit(self) -> None:
        """
        Clear fields related to previous run and change state to "ready".
        """
        status_from = self.data.status
        if status_from not in ("finished", "cancelled"):
            raise RunInputError(
                "Cannot resubmit run while previous run is still active."
            )

        logger.info(f"Run {self.data.id}: Resubmit run")
        timestamp = datetime.utcnow()
        log_entry = models.Run_Log(
            run_id=self.data.id,
            status_from=status_from,
            status_to="ready",
            timestamp=timestamp,
            reason="resubmit run",
        )
        try:
            savepoint = self.session.begin_nested()
            self.data.status = "ready"
            self.data.updated = timestamp
            self.data.workflow_root = None
            self.data.workflow_name = None
            self.data.cromwell_run_id = None
            self.data.result = None
            self.session.add(log_entry)
            self.session.commit()
        except SQLAlchemyError as error:
            savepoint.rollback()
            logger.exception(f"Unable to update Run {self.data.id}: {error}")
        return {self.data.id: "ready"}

    @retry(
        reraise=True,
        stop=stop_after_attempt(JAWS_GET_METADATA_MAX_RETRIALS),
        before=before_log(logger, logging.DEBUG),
        after=after_log(logger, logging.DEBUG),
        wait=wait_fixed(JAWS_GET_METADATA_WAIT_SEC),
    )
    def get_metadata(self, **kwargs):
        """
        Get Cromwell metadata, save for future use.
        Returns cached object unless "force" option is provided.
        A CromwellError exception may be raised if Cromwell is unreachable.
        """
        if self.data.cromwell_run_id is not None:
            force = kwargs.get("force", False)
            if force or self._metadata is None:
                try:
                    self._metadata = cromwell.get_metadata(self.data.cromwell_run_id)
                except Exception as e:
                    logger.critical(
                        f"cromwell.get_metadata raised an exception: {e} {self._metadata}"
                    )
                    raise CromwellGetMetadataError(
                        f"cromwell.get_metadata raised an exception: {e}"
                    )
            return self._metadata
        else:
            return None

    def check_cromwell_metadata(self):
        """
        Check Cromwell metadata for workflow_root and workflow_name.
        If no metadata found for a run id, return None
        :return: our JAWS Cromwell Metadata object
        :rtype: Cromwell.Metadata
        """
        metadata = None
        workflow_name = None
        workflow_root = None
        logger.debug(f"Run {self.data.id}: Check Cromwell Run metadata")
        try:
            metadata = self.get_metadata()
        except CromwellGetMetadataError as e:
            logger.warning(f"Can't find a metadata for {self.data.id}: {e}.")
            raise
        except Exception as e:
            logger.critical(f"Cromwell raises an exception {e}.")

        if metadata is not None:
            workflow_name = metadata.get("workflowName")
            workflow_root = metadata.get("workflowRoot")
        if workflow_name or workflow_root:
            try:
                self.data.workflow_name = workflow_name
                self.data.workflow_root = workflow_root
                self.session.commit()
            except SQLAlchemyError as error:
                self.session.rollback()
                logger.exception(f"Unable to update Run {self.data.id}: {error}")
        return metadata

    def check_cromwell_run_status(self) -> None:
        """
        Check Cromwell for the status of the Run.
        """
        if self.data.workflow_root is None:
            try:
                metadata = self.check_cromwell_metadata()
            except CromwellGetMetadataError as error:
                logger.error(
                    f"Run {self.data.id}: Failed to generate metadata: {error}"
                )
        try:
            cromwell_status = cromwell.get_status(self.data.cromwell_run_id)
        except CromwellError as error:
            logger.error(
                f"Run {self.data.id}: Unable to check Cromwell status: {error}"
            )
            return
        if not cromwell_status:
            return

        # check if state has changed.
        if cromwell_status == "Running":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.update_run_status("queued")
            elif self.data.status == "queued":
                # although Cromwell may consider a Run to be "Running", since it does not distinguish between
                # "queued" and "running", we check the task-log to see if any task is "running"; only once any
                # task is running does the Run transition to the "running" state.
                if self.did_run_start() is True:
                    self.update_run_status("running")
        elif cromwell_status == "Failed":
            err_msg = "Cromwell execution failed"
            if self.data.workflow_root is None:
                # The run failed in input processing stage, so no output folder was created, and the
                # errors report shall not be returned to the user.
                # Thus, we need to include the error message in the run-log.
                try:
                    metadata = self.get_metadata()
                except CromwellServiceError as error:
                    logger.error(
                        f"Run {self.data.id}: Failed to generate metadata: {error}"
                    )
                err_msg = "Cromwell submission failed"
                if metadata is not None:
                    errors_report = metadata.errors()
                    if (
                        "failures" in errors_report
                        and "message" in errors_report["failures"]
                    ):
                        err_msg = errors_report["failures"]
                    err_detail = []
                    if "causedBy" in errors_report:
                        for item in errors_report["causedBy"]:
                            err_detail.append(item["message"])
                    if len(err_detail):
                        err_msg = f"{err_msg}: " + "; ".join(err_detail)
            self.update_run_status(
                "failed", err_msg, cromwell_run_id=self.data.cromwell_run_id
            )
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.update_run_status("queued")
            if self.data.status == "queued":
                self.update_run_status("running")
            self.update_run_status(
                "succeeded", cromwell_run_id=self.data.cromwell_run_id
            )
        elif cromwell_status == "Aborted":
            self.update_run_status("cancelled")

    def update_run_status(self, status_to, reason=None, cromwell_run_id=None) -> None:
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        if reason is not None:
            reason = reason[:MAX_ERROR_STRING_LEN]
        status_from = self.data.status
        logger.info(f"Run {self.data.id}: now {status_to}")
        timestamp = datetime.utcnow()
        log_entry = models.Run_Log(
            run_id=self.data.id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
            cromwell_run_id=cromwell_run_id,
        )
        try:
            savepoint = self.session.begin_nested()
            self.data.status = status_to
            self.data.updated = timestamp
            if status_to in ("succeeded", "failed", "cancelled"):
                self.data.result = status_to
            self.session.add(log_entry)
            self.session.commit()
        except SQLAlchemyError as error:
            savepoint.rollback()
            logger.exception(f"Unable to update Run {self.data.id}: {error}")

    def write_supplement(self):
        """
        After Cromwell completes, several supplementary files are written to the Run's
        workflow_root folder.
        """
        logger.info(f"Run {self.data.id}: Write supplementary files")

        # get Cromwell metadata and set workflow_root if undefined
        try:
            metadata = self.check_cromwell_metadata()
        except CromwellGetMetadataError as e:
            logger.error(f"Run {self.data.id}: Failed to get metadata: {e}")
            return
        except CromwellServiceError as error:
            logger.error(f"Run {self.data.id}: Failed to generate metadata: {error}")
            self.update_run_status(
                "complete", "Cromwell metadata could not be retrieved"
            )
            return

        # confirm workflow_root is defined
        root = self.data.workflow_root
        if root is None:
            self.update_run_status("complete", "Cromwell run folder was not created")
            return

        # confirm the run folder exists
        if not os.path.isdir(root):
            self.update_run_status("complete", "Cromwell run folder does not exist")
            return

        # write metadata json
        metadata_file = os.path.join(root, "metadata.json")
        logger.debug(f"Run {self.data.id}: Writing {metadata_file}")
        write_json_file(metadata_file, metadata.data)

        # copy wdl
        infiles = []
        src_wdl_path = os.path.join(
            self.config["inputs_dir"], f"{self.data.submission_id}.wdl"
        )
        wdl_path = os.path.join(root, "workflow.wdl")
        if self.data.wdl_basename is not None:
            wdl_path = os.path.join(root, self.data.wdl_basename)
        try:
            shutil.copy(src_wdl_path, wdl_path)
        except IOError as error:
            logger.error(
                f"Run {self.data.id}: Failed to copy WDL to output dir: {error}"
            )
        else:
            infiles.append(self.data.wdl_basename)

        # copy subworkflows zip (if exists)
        src_subworkflows_path = os.path.join(
            self.config["inputs_dir"], f"{self.data.submission_id}.zip"
        )
        subworkflows_path = os.path.join(root, "subworkflows.zip")
        if os.path.isfile(src_subworkflows_path):
            try:
                shutil.copy(src_subworkflows_path, subworkflows_path)
            except IOError as error:
                logger.error(
                    f"Run {self.data.id}: Failed to copy subworkflows-ZIP to output dir: {error}"
                )
            else:
                infiles.append("subworkflows.zip")

        # copy inputs json
        src_inputs_json_path = os.path.join(
            self.config["inputs_dir"], f"{self.data.submission_id}.json"
        )
        inputs_json_path = os.path.join(root, "inputs.json")
        if self.data.json_basename is not None:
            inputs_json_path = os.path.join(root, self.data.json_basename)
        try:
            shutil.copy(src_inputs_json_path, inputs_json_path)
        except IOError as error:
            logger.error(
                f"Run {self.data.id}: Failed to copy inputs-JSON to output dir: {error}"
            )
        else:
            infiles.append(self.data.json_basename)

        # write errors report
        try:
            errors_report = metadata.errors()
        except Exception as error:
            logger.error(
                f"Run {self.data.id}: Failed to generate errors report: {error}"
            )
            self.update_run_status("complete", "Failed to generate errors report")
        else:
            errors_file = os.path.join(root, "errors.json")
            logger.debug(f"Run {self.data.id}: Writing {errors_file}")
            write_json_file(errors_file, errors_report)

        # write outputs
        try:
            outputs = metadata.outputs(relpaths=True)
        except Exception as error:
            logger.error(
                f"Run {self.data.id}: Failed to generate outputs file: {error}"
            )
            self.update_run_status("complete", "Failed to generate outputs file")
        else:
            outputs_file = os.path.join(root, "outputs.json")
            logger.debug(f"Run {self.data.id}: Writing {outputs_file}")
            write_json_file(outputs_file, outputs)

        # write outfiles
        try:
            outfiles = metadata.outfiles()
        except Exception as error:
            logger.error(
                f"Run {self.data.id}: Failed to generate outfiles file: {error}"
            )
            self.update_run_status("complete", "Failed to generate outfiles file")
        else:
            outfiles_file = os.path.join(root, "outfiles.json")
            logger.debug(f"Run {self.data.id}: Writing {outfiles_file}")
            write_json_file(outfiles_file, outfiles)

        # update task log
        try:
            task_log = self.task_log()
            task_summary_dict = metadata.task_summary_dict()
            task_log.add_metadata(task_summary_dict)
        except Exception as error:
            logger.error(f"Run {self.data.id}: Failed to update task-log: {error}")
            self.update_run_status("complete", "Failed to update task-log")

        # write task log
        try:
            task_log_table = task_log.table()
        except Exception as error:
            logger.error(f"Run {self.data.id}: Failed to generate task-log: {error}")
            self.update_run_status("complete", "Failed to generate task-log")
        else:
            task_log_file = os.path.join(root, "tasks.json")
            logger.debug(f"Run {self.data.id}: Writing {task_log_file}")
            write_json_file(task_log_file, task_log_table)

        # add cpu-hours to run summary
        self.data.cpu_hours = task_log.cpu_hours()

        # write run summary
        try:
            summary = self.summary()
        except Exception as error:
            logger.error(f"Run {self.data.id}: Failed to generate summary: {error}")
            self.update_run_status("complete", "Failed to generate summary")
        else:
            summary_file = os.path.join(root, "summary.json")
            logger.debug(f"Run {self.data.id}: Writing {summary_file}")
            write_json_file(summary_file, summary)

        # write output manifest (i.e. files to return to user)
        try:
            failed_folders = metadata.failed_folders()
        except Exception as error:
            logger.error(
                f"Run {self.data.id}: Failed to generate output manifest: {error}"
            )
            self.update_run_status("complete", "Failed to generate output manifest")
            return

        manifest = [
            *infiles,
            *outfiles,
            *failed_folders,
            "metadata.json",
            "errors.json",
            "outputs.json",
            "output_manifest.json",
            "summary.json",
            "tasks.json",
        ]
        manifest_file = os.path.join(root, "output_manifest.json")
        logger.debug(f"Run {self.data.id}: Writing {manifest_file}")
        write_json_file(manifest_file, manifest)
        self.update_run_status("complete")

    def output_manifest(self) -> list:
        """
        Return a list of the output files for a completed run, as paths relative to the workflow_root.
        """
        manifest_file = f"{self.data.workflow_root}/output_manifest.json"
        if not os.path.isfile(manifest_file):
            logger.error(f"Run {self.data.id}: Output manifest does not exist")
            return []
        try:
            manifest = read_json(manifest_file)
        except RunFileNotFoundError as error:
            logger.error(f"Run {self.data.id}: Failed to read output manifest: {error}")
            return []
        else:
            return manifest

    def publish_report(self):
        """
        Save final run metadata and send report document to reports service via RPC.
        We currently record resource metrics for successful and failed, but not cancelled Runs.
        """
        logger.info(f"Publish report for run {self.data.id}")

        if self.data.result == "cancelled":
            self.update_run_status("finished")
            return

        # read previously generated summary
        summary_file = f"{self.data.workflow_root}/tasks.json"
        if not os.path.isfile(summary_file):
            self.update_run_status("finished", "No report to publish to perf-metrics")
            return
        summary = self._read_json_file(summary_file)

        # publish and if successful, mark as done, else skip until next time
        try:
            response = self.reports_rpc_client.request(summary)
        except Exception as error:
            logger.exception(f"RPC save_run_report error: {error}")
            return
        if "error" in response:
            logger.warning(
                f"RPC save_run_report failed: {response['error']['message']}"
            )
            # do not change state; try again next time
        else:
            self.update_run_status("finished")


def check_active_runs(session, central_rpc_client, reports_rpc_client) -> None:
    """
    Get active runs from db and have each check and update their status.
    """
    active_states = [
        "ready",
        "submitted",
        "queued",
        "running",
        "succeeded",
        "failed",
        "cancel",
        "cancelled",
        "complete",
    ]
    try:
        rows = (
            session.query(models.Run).filter(models.Run.status.in_(active_states)).all()
        )
    except SQLAlchemyError as error:
        logger.warning(f"Failed to select active runs from db: {error}", exc_info=True)
    n = len(rows)
    if n == 0:
        return

    # If there is an unexpected error prevening a Run from being processed, it could block
    # other runs from being processed, so we randomize the list.  This isn't usually necessary,
    # but it increases the robustness of the system.
    if n > 1:
        shuffle(rows)
    for row in rows:
        run = Run(
            session,
            row,
            central_rpc_client=central_rpc_client,
            reports_rpc_client=reports_rpc_client,
        )
        run.check_status()


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
        run = Run.from_id(session, log.run_id)
        if run.data.cromwell_run_id is not None:
            data["cromwell_run_id"] = run.data.cromwell_run_id
        if run.data.workflow_root is not None:
            data["workflow_root"] = run.data.workflow_root
            data["workflow_name"] = run.data.workflow_name
        if log.status_to == "complete":
            data["output_manifest"] = run.output_manifest()
            data["cpu_hours"] = run.data.cpu_hours
        elif log.status_to in ("succeeded", "failed", "cancelled"):
            data["result"] = run.data.result
        try:
            response = central_rpc_client.request("update_run_log", data)
        except Exception as error:
            logger.exception(f"RPC update_run_log error: {error}")
            continue
        if "error" in response:
            logger.info(f"RPC update_run_status failed: {response['error']['message']}")
            continue
        try:
            log.sent = True
            session.commit()
        except Exception as error:
            session.rollback()
            logger.exception(f"Error updating run_logs as sent: {error}")


def max_active_runs_exceeded(
    session: sessionmaker, user_id: str, max_active_runs: int
) -> bool:
    """
    Get active runs from db and have each check and update their status.

    Returns True if max active runs exceeded, else False.
    """
    if max_active_runs is None or max_active_runs <= 0:
        return False
    active_states = [
        "submitted",
        "queued",
        "running",
    ]
    try:
        rows = (
            session.query(models.Run)
            .filter(models.Run.status.in_(active_states))
            .filter(models.Run.user_id == user_id)
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to get user active runs from db: {error}", exc_info=True
        )
        return True

    if len(rows) > max_active_runs:
        return True

    return False


class RunLog:
    """Class representing table of run state transition logs"""

    def __init__(self, session, run_id) -> None:
        self.session = session
        self.run_id = run_id
        self.data = self._select_rows()

    def _select_rows(self):
        try:
            rows = (
                self.session.query(models.Run_Log)
                .filter(models.Run_Log.run_id == self.run_id)
                .order_by(models.Run_Log.id)
                .all()
            )
        except SQLAlchemyError as error:
            logger.error(f"Unable to select from run_logs: {error}")
            raise RunDbError(error)
        return rows

    def logs_table(self):
        return self.data

    def logs(self):
        """Reformat logs table to (verbose) dictionary"""
        logs = []
        for row in self.data:
            (id, run_id, status_from, status_to, timestamp, reason, sent) = row
            log = {
                "status_from": status_from,
                "status_to": status_to,
                "timestamp": timestamp,
                "reason": reason,
            }
            logs.append(log)
        return logs


def s3_parse_uri(full_uri):
    """
    Extract the bucket name and object key from full URI
    :param full_uri: String containing bucket and obj key, starting with "s3://"
    :ptype full_uri: str
    :return: s3 bucket name, object key
    :rtype: list
    """
    full_uri = full_uri.replace("s3://", "", 1)
    folders = full_uri.split("/")
    s3_bucket = folders.pop(0)
    obj_key = "/".join(folders)
    return s3_bucket, obj_key


def read_json(path) -> any:
    with open(path, "r") as fh:
        contents = json.load(fh)
    return contents
