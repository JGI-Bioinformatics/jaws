"""
This class provides an API to a remotely running jaws-worker, via RPC.
"""

import logging
from datetime import datetime

# import sqlalchemy.exc
# from datetime import datetime
from jaws_site import models
from jaws_site.database import Session
from jaws_rpc.rpc_index import rpc_index


# config and logging must be initialized before importing this module
logger = logging.getLogger(__package__)


class WorkerDbError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WorkerGridError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WorkerRpcError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WorkerRpcOperationFailed(Exception):
    def __init__(self, message):
        super().__init__(message)


class WorkerNotFoundError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Worker(object):
    def __init__(self, worker_id=None):
        """
        Init Worker object.  If worker_id is specified, load from db.
        """
        self.worker_id = None
        if worker_id:
            self.worker_id = worker_id
            self.__retrieve()

    def __retrieve(self):
        """
        Load Worker values from db.
        """
        try:
            session = Session()
            worker = (
                session.query(models.Worker).filter_by(id=self.worker_id).one_or_none()
            )
        except Exception as error:
            raise WorkerDbError(f"Db unavailable: {error}")
        if not worker:
            raise WorkerNotFoundError
        self.worker = worker

    def __start(self, params):
        """
        Create a new Worker, submit to scheduler, and insert into db.
        """
        # validate input
        required_params = (
            "cpu",
            "max_time",
            "memory_gb",
            "job_id",
        )
        missing_params = []
        for required_param in required_params:
            if required_param not in params:
                missing_params.append(required_param)
        if missing_params:
            raise ValueError(f"Missing required fields: {', '.join(missing_params)}")

        # calculate max time in minutes (from string in "hh:mm:ss" format)
        max_time = params["max_time"].split(":")
        for n in range(len(max_time), 3):
            max_time.insert(0, 0)  # add any missing fields
        max_minutes = max_time[0] * 60 + max_time[1]

        # use the Site's datetime since the db, scheduler, and/or cluster nodes
        # could have wrong time (we require consistency for timedelta)
        created_timestamp = datetime.now()

        # insert into db, get worker_id (which is used to name job)
        session = Session()
        try:
            worker = models.Worker(
                job_id=params["job_id"],
                script=params["script"],
                job_name=params["job_name"],
                cwd=params["cwd"],
                out=params["out"],
                err=params["err"],
                max_time=params["max_time"],
                max_minutes=max_minutes,
                memory_gb=params["memory_gb"],
                status="created",
                pid=None,
                created_timestamp=created_timestamp,
            )
            session.add(worker)
            session.commit()
        except Exception as error:
            logger.exception(f"Failed to insert Worker: {error}")
            raise WorkerDbError(f"Db error: {error}")

        # save
        self.worker_id = worker.id
        self.worker = worker

    def stop(self):
        """Stop the worker."""
        raise NotImplementedError

    def minutes_remaining(self):
        """Return the time the worker has remaining before it stops."""
        start = self.worker.created_timestamp
        if start is None:
            return None
        delta = datetime.now() - start
        return int(delta.total_seconds() / 60)

    def __get_rpc_client(self):
        """
        Get RPC client from rpc_index
        """
        try:
            self.rpc = rpc_index.get_client(self.worker_id)
        except Exception as error:
            raise WorkerRpcError(f"Get RPC client failed: {error}")

    def __call(self, operation: str, params: dict) -> dict:
        """
        Send RPC call to worker.
        :param operation: Name of operation to request
        :type operation: str
        :param params: Parameters; varies by operation
        :type params: dict
        :return: JSON-RPC2 response; result content varies by operation
        :rtype: dict
        """
        try:
            response = self.rpc.request(operation, params)
        except Exception as error:
            logger.error(f"RPC {operation} failed: {error}")
            raise WorkerRpcError(f"RPC {operation} failed: {error}")
        return response

    def status(self) -> bool:
        """
        Check if jaws-worker is running and reachable via RPC.
        """
        response = self.__call("status")
        if "error" in response:
            return False
        return True

    def run(self, params):
        """
        Submit a task for execution.
        """
        raise NotImplementedError

    def kill(self, pid) -> None:
        """
        Instruct worker to kill Task.
        :param pid: Process id on worker node.
        :type pid: int
        :return: JSON-RPC2 response
        :rtype: dict
        """
        response = self.__call("kill", pid)
        if "error" in response:
            raise WorkerRpcOperationFailed(response["error"]["message"])

    def check_alive(self, pid):
        """
        Check if a task is running.
        """
        response = self.__call("check_alive", pid)
        if "error" in response:
            raise WorkerRpcOperationFailed(response["error"]["message"])
