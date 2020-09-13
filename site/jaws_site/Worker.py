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
            self.__load()

    def __load(self):
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

    def create(self, params: dict) -> int:
        """
        Create a new Worker, insert into db, and return worker_id.
        """
        # validate input
        required_params = (
            "job_id",
            "batch_name",
            "cwd",
            "max_minutes",
            "created_timestamp",
        )
        missing_params = []
        for required_param in required_params:
            if required_param not in params:
                missing_params.append(required_param)
        if missing_params:
            raise ValueError(f"Missing required fields: {', '.join(missing_params)}")

        # insert into db
        session = Session()
        try:
            worker = models.Worker(
                job_id=params["job_id"],
                batch_name=params["batch_name"],
                max_minutes=params["max_minutes"],
                status="queued",
                created_timestamp=params["created_timestamp"],
                array_task_id=None,  # will be added at start
                pid=None,  # will be added at start
            )
            session.add(worker)
            session.commit()
        except Exception as error:
            logger.exception(f"Failed to insert Worker: {error}")
            raise WorkerDbError(f"Db error: {error}")

        # save
        self.worker_id = worker.id
        self.worker = worker
        return self.worker_id

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
        try:
            self.__call("status")
        except Exception as error:
            logger.debug(f"Worker RPC status failed: {error}")
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
