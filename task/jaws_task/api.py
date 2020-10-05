"""
Classes associated with the JAWS Task Service.
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


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)

# constants
task_status_msg = {
    "ready": "The job has been prepared by Cromwell and submitted to JTM",
    "queued": "The job was received by JTM-manager and sent to JTM-worker",
    "pending": "The job was receive by JTM-worker and is awaiting resources",
    "running": "The job is currently executing",
    "success": "The job completed successfully",
    "failed": "The job has failed",
    "outofresource": "The job exceeeded the reserved RAM; increase the amount in the WDL and rerun",
    "terminated": "The run was terminated",
    "invalidtask": "The task definition in the message from jtm_submit is not valid",
    "timeout": "The worker timed out",
    "lostconnection": "The manager lost connection with the worker",
}


class DatabaseError(Exception):
    pass



## TODO CREATE TASK CLASS
class Task():

    def __init__(self, task_id, params=None):
        self.task_id = task_id

        # if new task, init new and save in db
        if params:
            self.save_task(params)

    def save_task():
        """
        Submit a task for execution.
        """
        # insert row into database
        pass  # TODO

    def status():
        """
        Get the status of a task.
        """
        # retrieve status from database
        pass  # TODO

    def cancel():
        """
        Abort the execution of a task.
        """
        # update row in database and get worker_id
        pass  # TODO

        # get rpc_client for worker
        pass  # TODO

        # send cancel instruction to jaws-worker via RPC
        pass # TODO


class TaskLog():

    def __init__(self):
        pass

def get_task_log(run_id):
    """Retrieve task log from database"""
    try:
        session = Session()
        query = (
            session.query(Job_Log).filter_by(run_id=run_id).order_by(Job_Log.timestamp)
        )
    except Exception as error:
        logger.exception(f"Error selecting task log: {error}")
        raise DatabaseError(f"{error}")

    result = []
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
        result.append(row)
    return result


def get_task_status(run_id):
    """
    Retrieve the current status of each task.
    """
    # get job log entries, sorted by timestamp,
    try:
        session = Session()
        query = (
            session.query(Job_Log).filter_by(run_id=run_id).order_by(Job_Log.timestamp)
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting from job_log: {error}")
        raise DatabaseError(f"{error}")

    # keep only the latest log entry per task
    tasks = collections.OrderedDict()
    for log in query:
        task_name = log.task_name
        reason = log.reason if log.reason else ""
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
    return result


def update_job_status(cromwell_run_id, cromwell_job_id, status_from, status_to, timestamp, reason):
    """
    A JTM worker shall post changes in job state, although it is missing the JAWS run id.
    """
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
        raise ValueError(f"{error}")

    # INSERT OR IGNORE because JTM sometimes sends duplicate messages
    session = Session()
    result = ""
    try:
        session.add(job_log)
        session.commit()
        logger.debug(f"Job {cromwell_job_id} status saved")
    except sqlalchemy.exc.IntegrityError:
        # duplicate message; ignore
        session.rollback()
    except Exception as error:
        session.rollback()
        session.close()
        raise DatabaseError(f"{error}")
    session.close()

    # if task is complete, transfer task output via Globus
    if status_to == "success":
        # TODO
        self.__transfer_folder(label, src_dir, dest_endpoint, dest_dir)

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
    rel_src_dir = os.path.relpath(src_dir, self.globus_default_dir)
    try:
        transfer_client = self.__authorize_transfer_client(transfer_rt)
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


class PriorityQueue():
    """
    A collection of tasks, retrievable by priority.
    Each worker pool has it's own priority queue.
    Every task in the queue is executable by any worker in the pool
    with sufficient time left to live.
    """

    def __init__(self, pool_id):
        pass

    def size(self):
        """
        Return the number of elements in the queue.
        """
        pass

    def add(self, task, run_id, factor):
        """
        Add a task to the queue with the parameters necessary to calculate it's priority.
        """
        pass

    def pop(self, worker_id, num_cpu, max_memory_gb, max_minutes, constraint=None):
        """
        Get the highest priority task which matches filter criteria, or None.
        """
        # TODO DON'T FORGET TO CAST RUN_ID AS FLOAT WHEN CALCULATING PRIORITY
        # SELECT MIN(CAST(run_id AS float) * factor)
        # may have to use two queries if using ORM
        pass

    def queue(self):
        """
        Get queued tasks, ordered by priority.
        """
        result = []
        return result

#    def stats(self):
#        """
#        Return number of tasks, total minutes requested, and maximum minutes requested.
#        """
#        num_tasks = 0
#        total_minutes = 0
#        max_minutes = 0
#        for task in self.tasks:
#            num_tasks = num_tasks + 1
#            total_minutes = total_minutes + task.max_minutes
#            if task.max_minutes > max_minutes:
#                max_minutes = task.max_minutes
#        return num_tasks, total_minutes, max_minutes


class WorkerPool():
    """
    A collection of workers and a priority queue.
    """

    def __init__(self, num_cpu, max_memory_gb, max_minutes, constraint=None):
        self.num_cpu = num_cpu
        self.max_memory_gb = max_memory_gb
        self.max_minutes = max_minutes
        self.constraint = constraint
        self.queue = PriorityQueue()
        self.workers = {}
        self.num_workers = 0

    def size(self):
        """
        Return the number of workers in the pool.
        """
        return self.num_workers

#    def stats(self):
#        """
#        Return number of workers, total minutes remaining, and maximum minutes remaining.
#        """
#        num_workers = 0
#        total_minutes = 0
#        max_minutes = 0
#        for worker in self.workers:
#            num_workers = num_workers + 1
#            total_minutes = total_minutes + worker.max_minutes
#            worker_max_minutes = worker.max_minutes()
#            if worker_max_minutes > max_minutes:
#                max_minutes = worker_max_minutes
#        return num_workers, total_minutes, max_minutes

    def evaluate(self):
        """
        Compare current queue to total workers.
        """
        queue = self.queue.queue()
        workers = self.workers()
        # run simulation to determine if worker pool size should be changed
        
         

    def increase(self, num_workers, max_minutes=480):
        """
        Increase the size of the worker pool by the specified amount.
        """
        pass

    def decrease(self, num_workers)
        """
        Decrease the size of the worker pool by the specified amount.
        """
        if num > self.num_workers:
            num = self.num_workers
        pass  # TODO


# TODO: This is just a db table to track active queues/worker pools
class Workforce():
    """
    A collection of Worker Pools.
    """

    def __init__(self):
        pass

    def add_task(self, task):
        """
        Add a task to the appropriate worker pool and return it.
        """
        pool = self.match(task)
        pool.add(task)
        return pool

    def filter(self, task):
        """
        Returns all worker pools which can run the task, sorted from smallest to largest.
        """

    def match(self, task):
        """
        Returns smallest worker pool which meets the requirements of the task.
        If necessary, the worker pool size is updated or a new worker pool is created.
        """
        for pool in self.pools.sorted():  # need custom sorted()
            if pool.match(task):
                pool.add(task)
                return pool

        # no matching pool was found; create anew
        pool = Pool(task)
        self.pools[pool.pool_id] = pool
        return pool

    def add(self, num_cpu, max_memory_gb, max_minutes, constraints={}):
        """
        Create a new worker pool and add it to the workforce.
        """
        pool = WorkerPool(num_cpu, max_memory_gb, max_minutes, constraint=None)
        constraintsList = []
        for key in constraints:
            constraintsList.append(f"{key}={constraints[key]}")
        constraints_sig = ';'.join(constraintsList.sorted())
        sig = f"{num_cpu}:{max_memory_gb}:{constraints_sig}"
        self.pools. 
        pass

    def retrieve(self, num_cpu, max_memory_gb, max_minutes, contraint=None):
        pass

    def update(self, delta):
        """
        """
        pass

    def delete(self, num_cpu, max_memory_gb, constraint=None):
        """
        Delete matching 
        """
        pass


class Daemon:
    """
    The daemon periodically checks priority-queue's and their worker-pools.
    """

    def __init__(self):
        """
        Init daemon with schedule
        """
        schedule.every(10).seconds.do(self.compare_queue_and_pool)

    def start_daemon(self):
        """
        Start the scheduled loop.
        """
        while True:
            schedule.run_pending()
            time.sleep(1)

    def compare_queue_and_pool(self):
        """
        Check on current status of active runs.  This is run as a new thread.
        """
        # since this is a new thread, must get new session
        session = db.Session()

        # get list of active runs
        try:
            active_queues = (
                db.session.query(db.Queue)
                .filter(db.Queue.status=="active")  # TODO
                .all()
            )
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting active queuess: {error}")
            raise DatabaseError(f"{error}")

        # init Queue objects and have them check their own worker-pools
        for row in active_queues:
            queue = queue.Queue(queue.id, session)
            queue.check_workers()

        session.close()
