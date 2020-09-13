"""
Task represents a task/job submitted by Cromwell.  It has a script and execution record.
"""

import logging
from jaws_site.database import Session
from jaws_site import models
from jaws_site.Worker import Worker


logger = logging.getLogger(__package__)


class Task(object):
    def __init__(self, task_id=None, params={}):
        """
        If task_id is provided, then load params from db.
        If params are provided, init new Task and save in db.
        """
        self.task_id = None
        if task_id and params:
            raise ValueError("Either task_id or params may be specified but not both.")
        elif task_id:
            self.task_id = task_id
            self.__retrieve()
        elif params:
            self.__create(params)

    def __retrieve(self):
        """
        Load Task values from db.
        """
        try:
            session = Session()
            task = session.query(models.Task).filter_by(id=self.task_id).one_or_none()
        except Exception as error:
            raise TaskDbError(f"Db unavailable: {error}")
        if not task:
            raise TaskNotFoundError
        self.task = task

    def __create(self, params):
        """
        Create a new Task and save in db.
        """
        # validate input
        required_params = (
            "script",
            "job_name",
            "cwd",
            "out",
            "err",
            "max_time",
            "memory_gb",
        )
        missing_params = []
        for required_param in required_params:
            if required_param not in params:
                missing_params.append(required_param)
        if missing_params:
            raise ValueError(f"Missing required fields: {', '.join(missing_params)}")

        # insert into db, get task_id
        try:
            session = Session()
            task = models.Task(
                script=params["script"],
                job_name=params["job_name"],
                cwd=params["cwd"],
                out=params["out"],
                err=params["err"],
                max_time=params["max_time"],
                memory_gb=params["memory_gb"],
            )
            session.add(task)
            session.commit()
        except Exception as error:
            logger.exception(f"Failed to insert Task: {error}")
            raise TaskDbError(f"Db error: {error}")

        # save in obj
        self.task_id = task.id
        self.task = task

    def check_alive(self) -> bool:
        """
        Returns True if alive, False otherwise.
        """
        if self.task.status in ("queued", "running"):
            return True
        else:
            return False

    def kill(self) -> None:
        """
        Cancel a Task.
        """
        if self.task.status == "queued":
            try:
                session = Session()
                self.task.status = "aborted"
                session.commit()
            except Exception as error:
                raise TaskDbError(f"Db error: {error}")
        elif self.task.status == "running":
            worker = Worker(self.task.worker_id)
            worker.kill(self.task.pid)


class TaskDbError(Exception):
    def __init__(self, message):
        super().__init__(message)


class TaskNotFoundError(Exception):
    def __init__(self, message):
        super().__init__(message)
