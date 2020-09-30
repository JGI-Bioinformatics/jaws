"""
jaws-task RPC operations
"""

import logging
import sqlalchemy.exc
from sqlalchemy.exc import SQLAlchemyError
import collections
from datetime import datetime
import os
import re
import globus_sdk
from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.db import Session, Job_Log
from jaws_rpc.responses import success, failure


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


def get_task_log(params):
    """Retrieve task log from database"""
    run_id = params["run_id"]
    try:
        session = Session()
        query = (
            session.query(Job_Log).filter_by(run_id=run_id).order_by(Job_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting from job_logs: {error}")
        return failure(500, f"Db error; {error}")
    except Exception as error:
        logger.exception(f"Error selecting from job_logs: {error}")
        return failure(500, f"Db error; {error}")
    table = []
    for log in query:
        # replace None with empty string
        reason = log.reason if log.reason else ""
        task_name = log.task_name
        if task_name is None:
            task_name = "<updating>"
        attempt = log.attempt
        if attempt is None:
            attempt = 1
        row = [
            task_name,
            attempt,
            log.cromwell_job_id,
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            reason,
        ]
        table.append(row)
    return success(table)


def get_task_status(params):
    """
    Retrieve the current status of each task.
    """
    run_id = params["run_id"]
    # get job log entries, sorted by timestamp,
    try:
        session = Session()
        query = (
            session.query(Job_Log).filter_by(run_id=run_id).order_by(Job_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting from job_log: {error}")
        return failure(500, f"Db error; {error}")
    tasks = collections.OrderedDict()
    for log in query:
        task_name = log.task_name
        reason = log.reason if log.reason else ""
        # keep only the latest log entry per task
        tasks[task_name] = [
            log.task_name,
            log.attempt,
            log.cromwell_job_id,
            log.status_from,
            log.status_to,
            log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
            reason,
        ]
        # this keeps the log entries sorted by timestamp
        tasks.move_to_end(task_name)

    # return only the values, the key (task_name) was just used for dereplication
    result = list(tasks.values())

    # add explanatory text column
    for row in result:
        status_to = row[4]
        row.append(jaws_constants.task_status_msg.get(status_to, ""))
    return success(result)



def update_job_status(params):
    """A JTM worker shall post changes in job state, although it is missing the JAWS run id.
    The state change is simply saved in the db; any other actions will be performed by the daemon."""
    cromwell_run_id = params["cromwell_run_id"]  # Cromwell's run/workflow UUID
    cromwell_job_id = params["cromwell_job_id"]  # JTM's task_id
    status_from = params["status_from"]
    status_to = params["status_to"]
    timestamp = datetime.strptime(params["timestamp"], "%Y-%m-%d %H:%M:%S")
    reason = None
    if "reason" in params:
        reason = params["reason"]
    if reason is not None:
        logger.info(
            f"Received job status: {cromwell_run_id}:{cromwell_job_id}:{status_from}:{status_to}:{reason}"
        )
    else:
        logger.info(
            f"Received job status: {cromwell_run_id}:{cromwell_job_id}:{status_from}:{status_to}"
        )

    # DEFINE ROW
    try:
        job_log = Job_Log(
            cromwell_run_id=cromwell_run_id,
            cromwell_job_id=cromwell_job_id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
        )
    except Exception as error:
        logger.exception(f"Failed to create job_log object for {params}: {error}")
        return failure(500, f"Failed to create job_log object for {params}: {error}")

    # INSERT OR IGNORE
    session = Session()
    result = ""
    try:
        session.add(job_log)
        session.commit()
        logger.debug(f"Job {cromwell_job_id} status saved")
    except sqlalchemy.exc.IntegrityError:
        # JTM sometimes sends duplicate messages; ignore
        logger.warning(f"Job {cromwell_job_id} status is duplicate; ignored")
        result = "Ignoring duplicate log entry"
    except Exception as error:
        session.rollback()
        session.close()
        logger.exception(f"Job {cromwell_job_id} status not saved: {error}")
        return failure(500, f"Failed to insert job {cromwell_job_id} log: {error}")
    session.close()
    return success(result)


def __authorize_transfer_client(self, token):
    client_id = config.conf.get("GLOBUS", "client_id")
    client = globus_sdk.NativeAppAuthClient(client_id)
    authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
    return globus_sdk.TransferClient(authorizer=authorizer)


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


def update_job_status_logs(self):
    """
    JTM job status logs are missing some fields: run_id, task_name, attempt.
    Fill them in now by querying Cromwell metadata.
    """
    logger.debug("Update job status logs")

    # select incomplete job log entries from database
    last_cromwell_run_id = None  # cache last
    last_metadata = None  # no need to get from cromwell repeatedly
    try:
        query = (
            self.session.query(Job_Log)
            .filter_by(task_name=None)
            .order_by(Job_Log.cromwell_run_id)
        )
    except Exception as error:
        logger.exception(f"Unable to select job_logs: {error}")
        return

    if not query:
        # nothing to do
        return

    for log in query:
        run_id = log.run_id  # may be None
        cromwell_run_id = log.cromwell_run_id
        cromwell_job_id = log.cromwell_job_id
        logger.debug(f"Job {cromwell_job_id} now {log.status_to}")

        # Lookup run_id given cromwell_run_id
        if not run_id:
            run = (
                self.session.query(Run)
                .filter_by(cromwell_run_id=cromwell_run_id)
                .one_or_none()
            )
            if not run:
                # JTM may send first job status update before cromwell_run_id is recorded
                continue
            log.run_id = run_id = run.id

            # update run state if first job to run
            if log.status_to == "running" and run.status == "queued":
                self.update_run_status(run, "running")

            self.session.commit()

        # TRY TO GET task_name AND attempt FROM CROMWELL METADATA
        if cromwell_run_id == last_cromwell_run_id:
            # another update for same run, reuse previously initialized object
            metadata = last_metadata
        else:
            # get from cromwell
            try:
                metadata = last_metadata = self.cromwell.get_metadata(cromwell_run_id)
                last_cromwell_run_id = cromwell_run_id
            except Exception as error:
                logger.exception(
                    f"Error getting metadata for {cromwell_run_id}: {error}"
                )
                continue
        job_info = metadata.get_job_info(cromwell_job_id)
        if job_info is None:
            status = metadata.get("status")
            if status in ("Failed", "Succeeded", "Aborted"):
                logger.error(
                    f"job_id {cromwell_job_id} not found in inactive run, {cromwell_run_id}"
                )
                # If the Run is done and the job_id cannot be found, it's an orphan job
                # (Cromwell never received the job_id from JTM), so set task_name to "ORPHAN"
                # to mark it as such and so it won't be picked up in the next round.
                log.task_name = "ORPHAN"
                self.session.commit()
                continue
            else:
                # The Cromwell metadata could be a bit outdated; try again next time
                logger.debug(
                    f"job_id {cromwell_job_id} not found in active run, {cromwell_run_id}"
                )
            continue
        log.attempt = job_info["attempt"]
        log.task_name = job_info["task_name"]
        task_dir = job_info[
            "call_root"
        ]  # not saved in db but may be used below for xfer
        self.session.commit()

        # if Task complete, then transfer output
        if log.status_to in ["success", "failed"] and task_dir:
            logger.debug(
                f"Transfer Run {log.run_id}, Task {log.task_name}:{log.attempt}"
            )

            # get Run record
            run = self.session.query(Run).get(log.run_id)

            # set task output dir
            task_output_dir = os.path.normpath(
                os.path.join(
                    run.output_dir,
                    os.path.relpath(task_dir, run.cromwell_workflow_dir),
                )
            )

            # set label
            short_task_name = log.task_name.split(".")
            short_task_name = short_task_name[1].split(":")
            label = f"Run {log.run_id} Task {short_task_name[0]}"
            label = re.sub("[^0-9a-zA-Z_]+", " ", label)

            # recursively transfer task dir
            transfer_task_id = self.__transfer_folder(
                label,
                run.transfer_refresh_token,
                task_dir,
                run.output_endpoint,
                task_output_dir,
            )
            log.debug(f"Xfer {log.task_name}: {transfer_task_id}")


# THIS DISPATCH TABLE IS USED BY jaws_rpc.rpc_server AND REFERENCES FUNCTIONS ABOVE
rpc_methods = {
    "update_job_status": {
        "method": update_job_status,
        "required_params": [
            "cromwell_run_id",
            "cromwell_job_id",
            "status_from",
            "status_to",
            "timestamp",
        ],
    },
    "get_task_log": {
        "method": get_task_log,
        "required_params": ["user_id", "run_id"],
    },
    "get_task_status": {
        "method": get_task_status,
        "required_params": ["user_id", "run_id"],
    },
}
