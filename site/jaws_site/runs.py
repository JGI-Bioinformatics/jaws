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


import os
import logging
import json
import io
from datetime import datetime
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.session import sessionmaker
import boto3
import botocore
from random import shuffle
from jaws_site import models, config
from jaws_site.cromwell import (
    Cromwell,
    CromwellError,
    CromwellRunNotFoundError,
    CromwellServiceError,
)
from jaws_site.tasks import TaskLog
from jaws_site.utils import write_json_file


logger = logging.getLogger(__package__)

cromwell = Cromwell(config.conf.get("CROMWELL", "url"))


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
        self.central_rpc_client = (
            kwargs["central_rpc_client"] if "central_rpc_client" in kwargs else None
        )
        self.reports_rpc_client = (
            kwargs["reports_rpc_client"] if "reports_rpc_client" in kwargs else None
        )

        try:
            self.config = {
                "site_id": config.conf.get("SITE", "id"),
                "inputs_dir": config.conf.get("SITE", "inputs_dir"),
                "default_container": config.conf.get(
                    "SITE", "default_container", "ubuntu:latest"
                ),
                "max_user_active_runs": int(
                    config.conf.get("SITE", "max_user_active_runs", 0)
                ),
                "aws_access_key_id": config.conf.get("AWS", "aws_access_key_id"),
                "aws_region_name": config.conf.get("AWS", "aws_region_name"),
                "aws_secret_access_key": config.conf.get(
                    "AWS", "aws_secret_access_key"
                ),
                "cromwell_url": config.conf.get("CROMWELL", "url"),
                "cromwell_executions_dir": config.conf.get(
                    "CROMWELL", "executions_dir"
                ),
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
            "workflow_root": self.data.workflow_root,
            "workflow_name": self.data.workflow_name,
            "submitted": self.data.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": self.data.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "status": self.data.status,
            "result": self.data.result,
            "compute_site_id": self.config["site_id"],
        }
        return summary

    def did_run_start(self):
        """
        Check if any task has started running by checking the task log.
        """
        if not self.data.cromwell_run_id:
            return False
        else:
            task_log = TaskLog(self.session, self.data.cromwell_run_id, logger)
            return task_log.did_run_start()

    def task_log(self):
        if not self.data.cromwell_run_id:
            return []
        task_log = TaskLog(self.session, self.data.cromwell_run_id)
        return task_log.table()

    def check_status(self) -> None:
        """Check the run's status, promote to next state if ready"""
        status = self.data.status
        logger.debug(f"Run {self.data.id} is {self.data.status}")
        if status in self.operations:
            return self.operations[status]()

    def mark_to_cancel(self) -> None:
        """
        Tag a run to be cancelled.  Raise if not successful.
        """
        if self.data.status in ["cancel", "cancelled"]:
            return
        try:
            self.update_run_status("cancel")
        except Exception as error:
            logger.error(f"Failed to mark Run {self.data.id} to be cancelled: {error}")
            raise RunDbError(
                f"Change Run {self.data.id} status to 'cancel' failed: {error}"
            )

    def cancel(self) -> None:
        """
        Cancel a run, aborting Cromwell if appropriate.

        Note: `cancel` = to-be-cancelled status
        """
        if self.data.cromwell_run_id and self.data.status in [
            "submitted",
            "queued",
            "running",
        ]:
            try:
                cromwell.abort(self.data.cromwell_run_id)
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
            raise IOError(error)
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
            raise IOError(f"File not found: {path}")
        data = None
        mode = "rb" if binary else "r"
        if not os.path.isfile(path):
            raise RunFileNotFoundError(f"File not found, {path}")
        try:
            with open(path, mode) as fh:
                data = fh.read()
        except IOError:
            raise
        if len(data) == 0:
            raise IOError("File is 0 bytes")
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
            return self._read_file_nfs(path, binary)

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
        with open(path, "w") as fh:
            fh.write(content)

    # TODO SWITCH FROM KWARGS TO CONF
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

    # TODO
    def read_inputs(self):
        """
        Read inputs json from file or S3 bucket and return contents.
        :return: The Run's input parameters
        :rtype: dict
        """
        fh = self._read_file(
            os.path.join(self.config["inputs_dir"], f"{self.data.submission_id}.json")
        )
        inputs = json.load(fh)
        return inputs

    def inputs(self):
        """
        Get the Run inputs with valid paths for this Site.
        :return: input parameters
        :rtype: dict
        """

        def add_prefix_to_paths(data, site_id, prefix):
            """Recursively traverse dictionary and add prefix to
            file path items"""
            if type(data) is str:
                if data.startswith(site_id):
                    return os.path.join(prefix, data)
                else:
                    return data
            elif type(data) is list:
                new_data = []
                for item in data:
                    new_data.append(add_prefix_to_paths(item, site_id, prefix))
                return new_data
            elif type(data) is dict:
                new_data = {}
                for key, value in data.items():
                    new_key = add_prefix_to_paths(key, site_id, prefix)
                    new_value = add_prefix_to_paths(value, site_id, prefix)
                    new_data[new_key] = new_value
                return new_data
            else:
                return data

        # convert paths to correct abspaths for this Site
        relpath_inputs = self.read_inputs()
        inputs = add_prefix_to_paths(
            relpath_inputs, self.data.input_site_id, self.config["inputs_dir"]
        )
        return inputs

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
        try:
            file_handles["inputs"] = self.inputs_fh()
        except Exception as error:
            raise RunFileNotFoundError(f"Error specifying inputs: {error}")
        try:
            path = os.path.join(
                self.config["inputs_dir"], f"{self.data.submission_id}.wdl"
            )
            file_handles["wdl"] = self._read_file(path)
        except Exception as error:
            raise RunFileNotFoundError(f"Cannot read {path}: {error}")
        try:
            path = os.path.join(
                self.config["inputs_dir"], f"{self.data.submission_id}.zip"
            )
            sub = self._read_file(path, True)
        except Exception:
            pass  # subworkflows are optional
        else:
            file_handles["subworkflows"] = sub
        return file_handles

    def cromwell_options(self):
        default_container = self.config["default_container"]
        options = {
            "caching": self.data.caching,
            "default_container": default_container,
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
        except Exception as error:
            self.update_run_status("submission failed", f"Input error: {error}")
        try:
            options = self.cromwell_options()
        except Exception as error:
            self.update_run_status("submission failed", f"Options error: {error}")
        try:
            cromwell_run_id = cromwell.submit(file_handles, options)
        except CromwellError as error:
            logger.error(f"Run {self.data.id} submission failed: {error}")
            # self.update_run_status("submission failed", f"{error}")
            return  # try again later
        else:
            self.data.cromwell_run_id = cromwell_run_id
            self.update_run_status("submitted")

    def resubmit(self) -> None:
        """
        Clear fields related to previous run and change state to "ready".
        """
        status_from = self.data.status
        if status_from != "finished":
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

    def check_cromwell_metadata(self) -> str:
        """
        Check Cromwell metadata for status, workflow_root, and workflow_name.
        This is used only for runs in "submitted" state.  Do not transition out of this
        state until we have the workflow_name and workflow_root.
        :return: status; may be None if unavailable
        :rtype: str
        """
        logger.debug(f"Run {self.data.id}: Check Cromwell Run metadata")
        try:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
        except CromwellError as error:
            logger.error(f"Run {self.data.id} failed to get Cromwell metadata: {error}")
            return None
        cromwell_status = metadata.get("status")
        workflow_name = metadata.get("workflowName")
        workflow_root = metadata.get("workflowRoot")
        if workflow_name is None or workflow_root is None:
            # setting these fields is a requirement to transition past this state
            return None
        logger.debug(
            f"Run {self.data.id} workflow_name={workflow_name}; workflow_root={workflow_root}"
        )
        try:
            self.data.workflow_name = workflow_name
            self.data.workflow_root = workflow_root
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Unable to update Run {self.data.id}: {error}")
            return None
        return cromwell_status

    def check_cromwell_run_status(self) -> None:
        """
        Check Cromwell for the status of the Run.
        """
        cromwell_status = None
        if self.data.status == "submitted":
            cromwell_status = self.check_cromwell_metadata()
        else:
            logger.debug(f"Run {self.data.id}: Check Cromwell Run status")
            try:
                cromwell_status = cromwell.get_status(self.data.cromwell_run_id)
            except CromwellError as error:
                logger.error(
                    f"Unable to check Cromwell status of Run {self.data.id}: {error}"
                )
                return
        if not cromwell_status:
            return

        # check if state has changed.
        logger.debug(f"Run {self.data.id}: Cromwell status is {cromwell_status}")
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
            self.update_run_status("failed")
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.update_run_status("queued")
            if self.data.status == "queued":
                self.update_run_status("running")
            self.update_run_status("succeeded")
        elif cromwell_status == "Aborted":
            self.update_run_status("cancelled")

    def update_run_status(self, status_to, reason=None) -> None:
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        status_from = self.data.status
        logger.info(f"Run {self.data.id}: now {status_to}")
        timestamp = datetime.utcnow()
        log_entry = models.Run_Log(
            run_id=self.data.id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
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
        # "test" is a special user account use to run end-to-end tests which do not require the
        # supplementary files, so skip.
        if self.data.user_id == "test":
            self.update_run_status("complete")
            self.update_run_status("finished")
            return

        logger.info(f"Run {self.data.id}: Write supplementary files")

        # get Cromwell metadata
        try:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
        except Exception as error:
            # TODO CHANGE EXCEPTIONS TO DISTINGUISH BETWEEN CONNECTION AND REJECTION ERRORS
            logger.warn(
                f"Run {self.data.id}: Unable to retrieve Cromwell metadata: {error}"
            )
            # TODO retry for x hrs before giving up and transitioning to finished
            self.update_run_status(
                "finished", "supplementary files could not be generated"
            )
            return

        # write metadata
        metadata_file = f"{self.data.workflow_root}/metadata.json"
        with open(metadata_file, "w") as fh:
            json.dump(metadata.data, fh)

        # write errors report
        errors_report = metadata.errors()
        errors_file = f"{self.data.workflow_root}/errors.json"
        write_json_file(errors_file, errors_report)

        # write outputs
        outputs = metadata.outputs(relpath=True)
        outputs_file = f"{self.data.workflow_root}/outputs.json"
        write_json_file(outputs_file, outputs)

        # write outfiles
        outfiles = metadata.outfiles(relpath=True)
        outfiles_file = f"{self.data.workflow_root}/outfiles.json"
        write_json_file(outfiles_file, outfiles)

        # write task summary
        task_summary = self.summary()
        # convert to old name so as not to upset Elasticsearch
        task_summary["site_id"] = task_summary["compute_site_id"].upper()
        del task_summary["compute_site_id"]
        task_summary["tasks"] = []
        for task in metadata.task_summary(last_attempts=True):
            # rename job_id key to cromwell_job_id
            task["cromwell_job_id"] = task["job_id"]
            task_summary["tasks"].append(task)
        summary_file = f"{self.data.workflow_root}/task_summary.json"
        write_json_file(summary_file, task_summary)

        # the task summary needs to be save in a db so it may be provided (e.g. to the dashboard)
        # even after the files have been purged
        try:
            self._insert_task_summary(task_summary)
        except Exception as error:
            logger.error(f"Run {self.data.id}: Unable to insert task-summary: {error}")
            # do not change state; try again later
            return

        self.update_run_status("complete")

    def _insert_task_summary(self, contents):
        task_summary = models.Task_Summary(
            run_id=self.data.id,
            tasks_json=contents,
        )
        try:
            savepoint = self.session.begin_nested()
            self.session.add(task_summary)
            self.session.commit()
        except SQLAlchemyError as error:
            savepoint.rollback()
            logger.exception(
                f"Unable to insert task-summary of Run {self.data.id}: {error}"
            )

    def publish_report(self):
        """
        Save final run metadata and send report document to reports service via RPC.
        We currently record resource metrics for successful and failed, but not cancelled Runs.
        """
        logger.info(f"Publish report for run {self.data.id}")

        if self.data.result == "cancelled":
            self.update_run_status("finished")
            return

        # read previously generate summary from file
        summary = self._read_json_file(f"{self.data.workflow_root}/tasks.json")

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
    logger.debug(f"There are {n} active runs")
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
        if log.status_to == "submitted":
            run = session.query(models.Run).get(log.run_id)
            data["cromwell_run_id"] = run.cromwell_run_id
        elif log.status_to == "queued":
            run = session.query(models.Run).get(log.run_id)
            data["workflow_root"] = run.workflow_root
            data["workflow_name"] = run.workflow_name
        try:
            response = central_rpc_client.request("update_run_logs", data)
        except Exception as error:
            logger.exception(f"RPC update_run_logs error: {error}")
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
