import os
import logging
import json
import io
from datetime import datetime
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm.exc import NoResultFound
import boto3
import botocore
from jaws_site import models, config
from jaws_site.cromwell import Cromwell, CromwellError

logger = logging.getLogger(__package__)

cromwell = Cromwell(config.conf.get("CROMWELL", "url"))


class RunDbError(Exception):
    pass


class RunNotFoundError(Exception):
    pass


class RunFileNotFoundError(Exception):
    pass


class Run:
    """Class representing a single Run"""

    def __init__(self, session, data, **kwargs):
        self.session = session
        self.data = data
        self.operations = {
            "ready": self.submit_run,
            "submitted": self.check_run_cromwell_status,
            "queued": self.check_run_cromwell_status,
            "running": self.check_run_cromwell_status,
            "succeeded": self.publish_report,
            "failed": self.publish_report,
        }
        self.central_rpc_client = (
            kwargs["central_rpc_client"] if "central_rpc_client" in kwargs else None
        )
        self.reports_rpc_client = (
            kwargs["reports_rpc_client"] if "reports_rpc_client" in kwargs else None
        )
        self._metadata = None

        self.config = {
            "site_id": config.conf.get("SITE", "id"),
            "inputs_dir": config.conf.get("SITE", "inputs_dir"),
            "default_container": config.conf.get(
                "SITE", "default_container", "ubuntu:latest"
            ),
            "aws_access_key_id": config.conf.get("AWS", "aws_access_key_id"),
            "aws_region_name": config.conf.get("AWS", "aws_region_name"),
            "aws_secret_access_key": config.conf.get("AWS", "aws_secret_access_key"),
            "cromwell_url": config.conf.get("CROMWELL", "url"),
        }

    @classmethod
    def from_params(
        cls, session, params, central_rpc_client=None, reports_rpc_client=None
    ):
        """Insert new Run into RDb.  Site only receives Runs in the "upload complete" state."""
        try:
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
            logger.warn(f"Run {run_id} not found: {error}")
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
            logger.warn(f"Cromwell {cromwell_run_id} not found: {error}")
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
        return self.data.status

    def summary(self) -> dict:
        """Produce summary of Run info"""
        summary = {
            "run_id": self.data.id,
            "user_id": self.data.user_id,
            "submitted": self.data.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": self.data.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "status": self.data.status,
            "result": self.data.result,
            "compute_site_id": self.config["site_id"],
        }
        return summary

    def metadata(self):
        """
        Lazy loading of Cromwell metadata.
        """
        if not self._metadata and self.data.cromwell_run_id:
            url = config.conf.get("CROMWELL", "url")
            self._metadata = Cromwell(url).get_metadata(self.data.cromwell_run_id)
        return self._metadata

    def did_run_start(self):
        return self.metadata().started_running()

    def task_summary(self):
        if self.data.cromwell_run_id:
            return self.metadata().task_summary()
        else:
            return []

    def task_log(self):
        if self.data.cromwell_run_id:
            return self.metadata().task_log()
        else:
            return []

    def report(self) -> dict:
        """Produce full report of Run and Task info"""
        # get cromwell metadata
        metadata = self.metadata()

        report = self.summary(last_attempt=True)
        report["run_id"] = self.data.id
        report["workflow_name"] = metadata.get("workflowName")
        report["cromwell_run_id"] = self.data.cromwell_run_id

        # self.summary returns a keyname compute_site_id. this was formely named site_id. To satisfy backwards
        # compatibility with kibana dashboard setup, need to rename compute_site_id back to site_id.
        report["site_id"] = report["compute_site_id"]
        del report["compute_site_id"]

        # transform task data structure and add to report
        report["tasks"] = []
        for task in self.task_summary():
            # delete some unwanted keys, rename others
            task["status"] = task["result"]
            del task["result"]
            del task["execution_status"]
            del task["job_id"]
            report["tasks"].append(task)

        return report

    def check_status(self) -> None:
        """Check the run's status, promote to next state if ready"""
        status = self.data.status
        if status in self.operations:
            return self.operations[status]()

    def cancel(self) -> None:
        """
        Cancel a run, aborting Cromwell or file transfer as appropriate.
        Raise if not successful.
        """
        orig_status = self.data.status
        try:
            self.update_run_status("cancelled")
        except Exception as error:
            logger.error(f"Failed to cancel Run {self.data.id}: {error}")
            raise RunDbError(
                f"Change Run {self.data.id} status to cancelled failed: {error}"
            )

        if self.data.cromwell_run_id and orig_status in [
            "submitted",
            "queued",
            "running",
        ]:
            try:
                cromwell.abort(self.data.cromwell_run_id)
            except CromwellError as error:
                logger.warn(f"Cromwell error cancelling Run {self.data.id}: {error}")
                raise CromwellError(f"Cromwell cancel failed: {error}")

    def outputs(self, relpath=True) -> dict:
        """
        Get outputs from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.data.cromwell_run_id:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            return metadata.outputs(relpath)
        else:
            return {}

    def outfiles(self, complete=False, relpath=True) -> dict:
        """
        Get output files from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.data.cromwell_run_id:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            return metadata.outfiles(complete=complete, relpath=relpath)
        else:
            return []

    def workflow_root(self) -> str:
        """
        Get workflowRoot from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.data.cromwell_run_id:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            return metadata.workflow_root()
        else:
            return None

    def output_manifest(self, complete=False) -> list:
        """
        Get list of all of a Run's output files.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.data.cromwell_run_id:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            result = {
                "workflow_root": metadata.workflow_root(),
                "manifest": metadata.outfiles(complete=complete, relpath=True),
            }
            return result
        else:
            return {}

    def errors(self) -> dict:
        """Get errors report from Cromwell service"""
        if self.data.cromwell_run_id:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
            return metadata.errors()
        else:
            return {}

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
        try:
            file_handles = self.get_run_inputs()
        except Exception as error:
            self.update_run_status("submission failed", f"Input error: {error}")
            return
        try:
            options = self.cromwell_options()
        except Exception as error:
            self.update_run_status("submission failed", f"Options error: {error}")
            return
        try:
            cromwell_run_id = cromwell.submit(file_handles, options)
        except CromwellError as error:
            logger.error(f"Run {self.data.id} submission failed: {error}")
            self.update_run_status("submission failed", f"{error}")
        else:
            self.data.cromwell_run_id = cromwell_run_id
            self.update_run_status("submitted", f"cromwell_run_id={cromwell_run_id}")

    def check_run_cromwell_status(self) -> None:
        """
        Check Cromwell for the status of the Run.
        """
        logger.debug(f"Run {self.data.id}: Check Cromwell status")
        try:
            cromwell_status = cromwell.get_status(self.data.cromwell_run_id)
        except CromwellError as error:
            logger.error(
                f"Unable to check Cromwell status of Run {self.data.id}: {error}"
            )
            raise

        # check if state has changed; allowed states and transitions for a successful run are:
        # submitted -> queued -> running -> succeeded (with no skipping of states allowed)
        # Additionally, any state can transition directly to cancelled or failed.
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
        self._update_run_status(status_to, timestamp)
        self._insert_run_log(status_from, status_to, timestamp, reason)

    def _update_run_status(self, new_status, timestamp) -> None:
        """
        Update Run's current status in 'runs' table
        """
        try:
            self.data.status = new_status
            self.data.updated = timestamp
            if new_status in ("succeeded", "failed"):
                self.data.result = new_status
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

    def publish_report(self):
        """
        Send report document to reports service via RPC.
        """
        # "test" is a special user account for automatic periodic system tests -- skip
        if self.data.user_id != "test":
            report = self.report()
            if not report:
                logger.exception("RPC save_run_report warning: run summary is empty.")
                return
            try:
                response = self.reports_rpc_client.request(report)
            except Exception as error:
                logger.exception(f"RPC save_run_report error: {error}")
                return
            if "error" in response:
                logger.warn(
                    f"RPC save_run_report failed: {response['error']['message']}"
                )
                return
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
    ]
    try:
        rows = (
            session.query(models.Run).filter(models.Run.status.in_(active_states)).all()
        )
    except SQLAlchemyError as error:
        logger.warning(f"Failed to select active runs from db: {error}", exc_info=True)
    else:
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
            (run_id, status_from, status_to, timestamp, reason, sent) = row
            log = {
                "status_from": status_from,
                "status_to": status_to,
                "timestamp": timestamp,
                "reason": reason,
            }
            logs.append(log)
        return logs
