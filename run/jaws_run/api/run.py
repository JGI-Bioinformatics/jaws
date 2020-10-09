"""
JAWS Analysis Service API

This service stores persistent Run information in a db and interacts with Cromwell.
"""

import logging
import os
import requests
from datetime import datetime
import sqlalchemy.exc
from sqlalchemy.exc import SQLAlchemyError
from jaws_run import config, db
from jaws_run.cromwell import Cromwell
from jaws_rpc.rpc_client import RpcClient, RpcError
from jaws_run.api.transfer import Transfer


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    pass


class RunNotFoundError(Exception):
    pass


class RunAlreadyExistsError(Exception):
    pass


class CromwellError(Exception):
    pass


class TaskServiceError(Exception):
    pass


class Run:
    """
    Representation of a run, which is the execution of a workflow (WDL) on a specific set of inputs (JSON) by Cromwell.
    A run's persistent data is stored in a relational database.
    If the Run does not already exist, the parameters to initialize it are required.
    """

    # init class vars
    client_id = config.conf.get("GLOBUS", "client_id")
    site_id = config.conf.get("SITE", "id")
    globus_root_dir = config.conf.get("GLOBUS", "root_dir")
    globus_default_dir = config.conf.get("GLOBUS", "default_dir")
    uploads_dir = os.path.join(
        config.conf.get("GLOBUS", "root_dir"),
        config.conf.get("SITE", "uploads_subdirectory"),
    )
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    downloads_dir = os.path.join(
        config.conf.get("GLOBUS", "root_dir"),
        config.conf.get("SITE", "downloads_subdirectory"),
    )
    globus_endpoint = config.conf.get("GLOBUS", "endpoint_id")
    downloads_subdir = config.conf.get("SITE", "downloads_subdirectory")
    terminal_states = ["cancelled", "download complete"]

    def __init__(self, run_id: int, session=None, params=None):
        """
        Initialize new object.  Since python doesn't support multiple constructors,
        either run_id is required (for an existing Run) or the required parameters
        (for a new Run) are required.
        """

        # init obj vars
        self.run_id = run_id
        self.data = None  # row in Runs table
        self.dispatch = {
            "created": self.__begin_upload,
            "uploading": self.__check_upload,
            "upload complete": self.__submit,
            "submitted": self.__check_cromwell_status,
            "queued": self.__check_cromwell_status,
            "running": self.__check_cromwell_status,
            "succeeded": self.__begin_download,
            "failed": self.__begin_download,
            "downloading": self.__check_download,
        }

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
            raise RunNotFoundError(f"{error}")
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting on runs table: {error}")
            raise DatabaseError(f"{error}")
        if not run:
            logger.debug(f"Run {self.run_id} not found")
            raise RunNotFoundError(f"No such record")
        self.data = run

    def get_status(self):
        """
        Return the current status of the Run.
        """
        return self.data.status

    def get_log(self):
        """
        Get complete Run Log.
        """
        try:
            logs = self.session.query(db.Run_Log).filter_by(run_id=self.run_id).all()
        except sqlalchemy.exc.IntegrityError as error:
            logger.exception(f"Run Logs not found for Run {self.run_id}: {error}")
            raise RunNotFoundError(f"{error}")
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting on run_logs table: {error}")
            raise DatabaseError(f"{error}")
        if not logs:
            logger.debug(f"Run {self.run_id} not found in logs")
            raise RunNotFoundError(f"No such record")
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
            logger.debug(
                f"Run {self.run_id} is {status}: Instructing Cromwell to cancel"
            )
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

    # methods which call jaws-task service:

    def get_task_service_client(self):
        if not self.task_service:
            self.task_service = RpcClient(config.conf.task_service_params())
        return self.task_service

    def get_task_log(self):
        """
        Retrieve log of all tasks' state transitions for a run.

        :return: table of task_id, task_name, attempt, status_from, status_to, timestamp, reason
        :rtype: list
        """
        logger.info(f"Get task log for Run {self.run_id}")
        if self.data.cromwell_run_id is None:
            return None
        try:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
        except requests.exceptions.HTTPError as error:  # TODO change this to CromwellConnectionError
            logger.exception(f"Unable to get metadata for {self.run_id}: {error}")
            raise CromwellError(f"{error}")

        # get tasks associated with this run
        tasks = metadata.tasks()  # TODO new method

        # make list of tasks
        task_ids = []
        tasks = {}  # task_id => Task object
        for task in metadata.tasks():
            # NB: cromwell's "job_id" = jaws' "task_id"
            task_ids.append(task.cromwell_job_id)
            tasks[task.cromwell_job_id] = task

        # have jaws-task service provide logs for requested tasks
        task_service = self.get_task_service_client()
        try:
            response = task_service("get_logs", {"task_ids": task_ids})
        except RpcClient.RpcError as error:
            raise RpcError(f"{error}")
        if "error" in response:
            raise TaskServiceError(f"Task service error: {response['error']['message']}")
        logs = response["result"]

        # combine cromwell task metadata and task log into table
        table = []
        for task_id in task_ids.sorted():
            task = tasks[task_id]
            for log in logs[task_id]:
                row = [task_id, task.name, task.attempt]
                row.append(log)
                table.append(row)
        return table

    def get_task_status(self):
        """
        Retrieve current status of each task of a run.

        :return: Table of task_id, task_name, attempt, status, timestamp, reason
        :rtype: list
        """
        logger.info(f"Get task status for Run {self.run_id}")
        if self.data.cromwell_run_id is None:
            return None
        try:
            metadata = cromwell.get_metadata(self.data.cromwell_run_id)
        except requests.exceptions.HTTPError as error:  # TODO change this to CromwellConnectionError
            logger.exception(f"Unable to get metadata for {self.run_id}: {error}")
            raise CromwellError(f"{error}")

        all_task_ids = []
        check_task_ids = {}  # tasks to request status from jaws-task service
        status = {}  # task_id => row
        for task in metadata.tasks():
            task_id = task.cromwell_job_id
            all_task_ids.append(task_id)

            if task.execution_status == "Done":
                done_status = "succeeded"
                reason = ""
                if task.return_code != 0:
                    done_status = "failed"
                    reason = task.error(task.attempt)
                timestamp = datetime.strptime(task.end, "%Y-%m-%dT%H:%M:%S.%fZ")
                status[task_id] = [
                    task_id,
                    task.name,
                    task.attempt,
                    done_status,
                    timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                    reason,
                ]
            else:
                # for these tasks, we shall get status from the Task service instead because
                # it has more detail than Cromwell
                status[task_id] = [task_id, task.name, task.attempt]
                check_task_ids[task_id] = None

        # get status from jaws-task service
        task_service = self.get_task_service_client()
        try:
            response = task_service("get_status", {"task_ids": check_task_ids.keys()})
        except RpcClient.RpcError as error:
            raise RpcError(f"{error}")
        if "error" in response:
            raise TaskServiceError(f"Task service error: {response['error']['message']}")
        check_status = response["result"]

        # combine cromwell task metadata and task log into table
        table = []
        for task_id in all_task_ids.sorted():
            task = status[task_id]
            row = status[task_id]
            if task_id in check_task_ids:
                row.append(check_status[task_id])
            table.append(row)
        return table

    def check_status(self):
        """
        Check if run is ready to transition to the next state.
        This method is called by the Daemon periodically and is the only way a Run
        transitions to a new state.
        """
#            "created": self.__begin_upload,
#            "uploading": self.__check_upload,
#            "upload complete": self.__submit,
#            "submitted": self.__check_cromwell_status,
#            "queued": self.__check_cromwell_status,
#            "running": self.__check_cromwell_status,
#            "succeeded": self.__begin_download,
#            "failed": self.__begin_download,
#            "downloading": self.__check_download,
        if self.status in self.dispatch:
            self.dispatch.get(self.status)()

    def __begin_upload(self):
        """
        A Run in "created" state shall transition to "uploading" upon submitting it's
        Globus transfer manifest to the Globus transfer service.  If unsuccessful, the Run
        transitions to "upload failed."
        """
        pass  # TODO

    def __get_upload(self):
        """
        Get upload Transfer object.
        """
        if self.data.upload_id:
            self.upload = Transfer(self.data.upload_id, self)
        else:
            self.upload = None
        return self.upload

    def __get_download(self):
        """
        Get download Transfer object.
        """
        if self.data.download_id:
            self.download = Transfer(self.data.download_id, self)
        else:
            self.download = None
        return self.download

    def __check_upload(self):
        """
        Check if Globus upload task is complete.
        A Run in "uploading" will transition to either "upload complete" or "upload failed".
        """
        logger.debug(f"Run {self.run_id}: Check upload status")
        upload_task_id = self.__get_upload()
        try:
            globus_status = self._get_globus_transfer_status(upload_task_id)
        except Exception as error:
            logger.exception(f"Failed to check upload {self.run_id}: {error}")
            raise
        if globus_status == "FAILED":
            self.__update_status__("upload failed")
        elif globus_status == "INACTIVE":
            self.__update_status__(
                "upload inactive",
                "Either your Globus endpoint is unavailable or your endpoint authorization has expired",
            )
        elif globus_status == "SUCCEEDED":
            self.__update_status__("upload complete")

    def __submit(self, run):
        """
        Submit the run to Cromwell.
        A Run in "upload complete" will transition to "submitted" if accepted by Cromwell service; or
        "submission failed" otherwise.
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
            self.__update_status__("missing input", "Missing inputs JSON file")
            return
        if not os.path.exists(wdl_file):
            logger.warning(f"Missing WDL for run {self.run_id}: {wdl_file}")
            self.__update_status__("missing input", "Missing WDL file")
            return
        if not os.path.exists(zip_file):
            zip_file = None

        # submit to Cromwell
        try:
            cromwell_run_id = cromwell.submit(wdl_file, json_file, zip_file)
        except Exception:
            raise
        if cromwell_run_id:
            self.data.cromwell_run_id = cromwell_run_id
            self.__update_status__(
                "submitted", f"cromwell_run_id={self.data.cromwell_run_id}"
            )
        else:
            self.__update_status__("submission failed")

    def __check_cromwell_status(self, run):
        """
        Check Cromwell for the status of one Run.
        A Run in "submitted" state will transition to "queued" once Cromwell reports it as
        "Running" because Cromwell doesn't have a queued state.
        In order to move a Run from "queued" to "running", the jaws-task service is queried.
        Ultimately, Cromwell will either report the status as "Succeeded" or "Failed", which
        shall move the Run to "succeeded" or "failed" state.
        It's also possible a Run could be cancelled, usually by the user.
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
        # Note: Cromwell and JAWS states are distinct and Cromwell lacks a "queued" state.
        logger.debug(f"Run {self.run_id}: Cromwell status is {cromwell_status}")
        if cromwell_status == "Running":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.__update_status__("queued")
        elif cromwell_status == "Failed":
            self.__update_status__("failed")
        elif cromwell_status == "Succeeded":
            # no skips allowed, so there may be more than one transition
            if self.data.status == "submitted":
                self.__update_status__("queued")
            if self.data.status == "queued":
                self.__update_status__("running")
            if self.data.status == "running":
                self.__update_status__("succeeded")
        elif cromwell_status == "Aborted":
            self.__update_status__("cancelled")

    def __begin_download(self, run, metadata=None):
        """
        Whether or not a Run ends up in "succeeded" or "failed" state, the Cromwell output is
        returned to the user.  The output was returned as each task completes, but we send the
        entire output again to be sure, using the Globus transfer option not to resend existing
        outfiles.
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
                self.__update_status__("download complete", "No run folder was created")
                return
        # ECCE
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
            self.__update_status__(
                "downloading", f"download_task_id={self.data.download_task_id}"
            )
            self.session.commit()
        # /ECCE
        # else ignore error; was logged; will try again later

    def __check_download(self, run):
        """
        If download is complete, change state to "download complete"; if transfer failed,
        then state becomes "download failed."
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
            self.__update_status__("download complete")
        elif globus_status == "FAILED":
            self.__update_status__("download failed")

    def __update_status__(self, new_status, reason=None):
        """
        Whenever a Run enters a new state, the transition is recorded in the run_log and
        the runs table is updated.
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
