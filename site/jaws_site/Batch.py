"""
A Batch is a collection of similar Workers (i.e. same cpu, ram).
Workers in a group are interchangeable; any one can perform a task assigned to the Batch,
granted they have sufficient time remaining.
"""

import logging
import datetime
from jaws_site import config

# from jaws_site import models
# from jaws_site.database import Session
from jaws_site.Worker import Worker

# from jaws_site.Task import Task
from jaws_site.GridFactory import GridFactory


logger = logging.getLogger(__package__)


class BatchGridError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Batch(object):
    def __init__(self, cpu, memory_gb, params):
        """
        Init object.
        """
        self.workers = []
        self.cpu = cpu
        self.memory_gb = memory_gb
        self.params = params  # optional params, e.g. account, qos

    def __load(self):
        """
        Read from db.
        """

    def status(self):
        """
        Summarize status of workers in batch.
        """
        result = {}
        for worker in self.workers:
            status = worker.status
            if status in result:
                result[status] = result[status] + 1
            else:
                result[status] = 1
        return result

    def evaluate_size(self):
        """
        Compare queue to current batch, adjust size as appropriate.
        """
        raise NotImplementedError

    def increase_size(self, num: int) -> None:
        """
        Add one or more workers matching the specified requirements.
        """
        # calculate max time in minutes (from string in "hh:mm:ss" format)
        max_time = self.params["max_time"].split(":")
        for n in range(len(max_time), 3):
            max_time.insert(0, 0)  # add any missing fields
        max_minutes = max_time[0] * 60 + max_time[1]

        # use the Site's datetime since the db, scheduler, and/or cluster nodes
        # could have wrong time (we require consistency for timedelta)
        created_timestamp = datetime.now()

        # submit to scheduler, get worker_ids
        params = {
            "array_size": num,  # will be ignore if not applicable
        }
        try:
            grid = GridFactory(config.conf("SITE", "scheduler"))
        except Exception as error:
            raise BatchGridError(f"{error}")
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
            job_ids = grid.submit(params)
        except Exception as error:
            raise BatchGridError(f"{error}")

        # initialize Worker objects, which are persistent
        # (i.e. stored in db)
        for job_id in job_ids:
            # job_id is str (e.g. "9999:1")
            worker = Worker()
            worker.create(job_id, created_timestamp, max_minutes)
            self.workers.append(worker)

    def decrease_size(self, num: int) -> None:
        """
        Release one or more workers, starting with oldest.
        """
        num_workers = len(self.workers)
        if num_workers - num < 0:
            num = num_workers
        for _ in range(num):
            worker = self.workers.pop(0)
            worker.stop()

    def stop(self) -> None:
        """
        Stop all workers.
        """
        for worker in self.workers:
            worker.stop()
        del self.workers[:]
