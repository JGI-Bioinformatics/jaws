"""
JAWS Analysis Service API

This service stores persistent Run information in a db and interacts with Cromwell.
"""

import logging
import os
import requests
import schedule
import time
from datetime import datetime
import sqlalchemy.exc
from sqlalchemy.exc import SQLAlchemyError
import globus_sdk
from jaws_run import config, db
from jaws_run.cromwell import Cromwell
from jaws_rpc.rpc_client import RpcClient, RpcError


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    """Generic exception when Cromwell does not return metadata."""

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
    No other object has direct access to Cromwell.
    A run's persistent data is stored in a relational database.
    If the Run does not already exist, the parameters to initialize it are required.
    """

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
        self.run_id = run_id
        self.data = None  # row in Runs table
        self.dispatch = {
            "created": self.begin_upload,
            "uploading": self.check_upload,
            "upload complete": self.submit,
            "submitted": self.check_cromwell_status,
            "queued": self.check_cromwell_status,
            "running": self.check_cromwell_status,
            "succeeded": self.begin_download,
            "failed": self.begin_download,
            "downloading": self.check_download,
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

    def get_upload(self):
        """
        Get upload Transfer object.
        """
        if self.data.upload_id:
            self.upload = Transfer(self.data.upload_id, self)
        else:
            self.upload = None
        return self.upload

    def get_download(self):
        """
        Get download Transfer object.
        """
        if self.data.download_id:
            self.download = Transfer(self.data.download_id, self)
        else:
            self.download = None
        return self.download

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
            response = task_service("get_statuses", {"task_ids": check_task_ids.keys()})
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

# TODO ECCE

    def check_status(self):
        """
        Check if run is ready to transition to the next state.
        """
        if self.status in self.dispatch:
            self.dispatch.get(self.status)()

    def begin_upload(self):
        """
        """
        pass  # TODO

    def check_upload(self):
        """
        Check if Globus upload task is complete.
        """
        logger.debug(f"Run {self.run_id}: Check upload status")
        upload_task_id = self.get_upload_task_id()
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

    def submit(self, run):
        """
        Submit the run to Cromwell.
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

    def check_cromwell_status(self, run):
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

    def begin_download(self, run, metadata=None):
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

    def check_download(self, run):
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
            self.__update_status__("download complete")
        elif globus_status == "FAILED":
            self.__update_status__("download failed")

    def __update_status__(self, new_status, reason=None):
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


class Transfer:
    """
    Class representing a Globus transfer task.
    """

    user_service_client = None  # rpc client for jaws-user service
    transfer_client = None  # globus transfer client obj

    def __init__(self, **kwargs):
        if "transfer_id" in kwargs.items:
            # an existing Globus transfer task
            self.transfer_id = kwargs.get("transfer_id")
            self.__load__()
        else:
            # a new Globus transfer task
            self.submit(**kwargs)

    def __load__(self):
        pass  # TODO

    def submit(self, **kwargs):
        """
        Submit a transfer to Globus, get transfer task id.
        """
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
        self.transfer_id = None  # TODO

    def user_service(self):
        if not self.user_service_client:
            self.user_service_client = RpcClient.client(
                config.conf.user_service_params()
            )
        return self.user_service_client

    def transfer_client(self):
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

    def status(self):
        """
        Query Globus transfer service for transfer task status.
        """
        transfer_client = self.get_transfer_client()
        try:
            task = transfer_client.get_task(self.transfer_id)
            self.status = task["status"]
        except Exception:
            logger.exception("Failed to check Globus upload status", exc_info=True)
            raise
        return self.status

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


class Daemon:
    """
    The daemon periodically checks on active runs, which may usher them to the next state.
    """

    active_run_states = [
        "uploading",
        "upload complete",
        "submitted",
        "queued",
        "running",
        "succeeded",
        "failed",
        "downloading",
    ]

    def __init__(self):
        """
        Init daemon with schedule
        """
        schedule.every(10).seconds.do(self.check_active_runs)

    def start_daemon(self):
        """
        Start the scheduled loop.
        """
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_active_runs(self):
        """
        Check on current status of active runs.  This is run as a new thread.
        """
        # since this is a new thread, must get new session
        session = db.Session()

        # get list of active runs
        try:
            active_runs = (
                db.session.query(db.Run)
                .filter(db.Run.status.in_(self.active_run_states))
                .all()
            )
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting active runs: {error}")
            raise DatabaseError(f"{error}")

        # init Run objects and have them check and update their status
        for row in active_runs:
            run = Run(row.id, session)
            run.check_status()

        session.close()
