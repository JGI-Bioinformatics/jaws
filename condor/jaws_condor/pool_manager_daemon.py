"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

from jaws_condor.pool_manager import PoolManager
import schedule
import time
import logging
from jaws_condor import config

logger = logging.getLogger(__package__)


class PoolManagerDaemon:
    """
    Daemon that periodically checks on this site's active runs.
    It prompts active runs to query Cromwell or Globus, as appropriate.
    When Globus uploads are completed, the run is submitted to Cromwell.
    When Cromwell execution has completed, the output is downloaded using Globus.
    When the download is finished, the run is complete.
    """

    def __init__(self):
        logger.info("Initializing pool manager daemon")
        self.time_add_worker_pool = config.conf.get("POOL_MANAGER", "time_add_worker_pool")
        self.time_rm_worker_pool = config.conf.get("POOL_MANAGER", "time_rm_worker_pool")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(self.time_add_worker_pool).minutes.do(self.add_worker_pool)
        schedule.every(self.time_rm_worker_pool).minutes.do(self.rm_worker_pool)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def add_worker_pool(self):
        """
        Check for runs in particular states.
        """
        pool = PoolManager()
        pool.add_workers()

    def rm_worker_pool(self):
        """
        Check for runs in particular states.
        """
        pool = PoolManager()
        pool.rm_workers()
