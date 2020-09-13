"""
This class provides an API to a remotely running jaws-worker, via RPC.
"""

import logging

# import sqlalchemy.exc
# from datetime import datetime
from jaws_site import config, models
from jaws_site.database import Session
from jaws_rpc.rpc_index import rpc_index
from jaws_grid.GridFactory import GridFactory


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
        )
        missing_params = []
        for required_param in required_params:
            if required_param not in params:
                missing_params.append(required_param)
        if missing_params:
            raise ValueError(f"Missing required fields: {', '.join(missing_params)}")

        # insert into db, get worker_id (which is used to name job)
        session = Session()
        try:
            worker = models.Worker(
                script=params["script"],
                job_name=params["job_name"],
                cwd=params["cwd"],
                out=params["out"],
                err=params["err"],
                max_time=params["max_time"],
                memory_gb=params["memory_gb"],
                status="created",
                job_id=None,
                pid=None,
            )
            session.add(worker)
            session.commit()
        except Exception as error:
            logger.exception(f"Failed to insert Worker: {error}")
            raise WorkerDbError(f"Db error: {error}")

        # save in obj
        self.worker_id = worker.id
        self.worker = worker

        # submit to scheduler
        try:
            scheduler = GridFactory(config.conf("SITE", "scheduler"))
        except Exception as error:
            raise WorkerGridError(f"{error}")
        params = {}
        #            "script":
        #            "job_name":
        #            "cwd":
        #            "out":
        #            "err":
        #            "qos":
        #            "time":
        #            "cpu":
        #            "constraint":
        #            "account":
        #            "exclusive":
        #            "memory":
        #            "mem-per-cpu":
        #        }
        try:
            self.job_id = scheduler.submit(params)
        except Exception as error:
            raise WorkerGridError(f"{error}")

        # save in db
        self.worker.job_id = self.job_id
        session.commit()

    def stop(self):
        """Stop the worker."""
        raise NotImplementedError

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
