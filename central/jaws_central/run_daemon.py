import schedule
import time
import logging
from jaws_central import database
from jaws_central.runs import check_active_runs


logger = logging.getLogger(__package__)


class RunDaemon:
    """
    Daemon that periodically checks on active runs.
    """

    def __init__(self, rpc_index):
        logger.info("Initializing run daemon")
        self.rpc_index = rpc_index

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.check_active_runs)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_active_runs(self):
        """
        Check for runs in particular states.
        """
        session = database.Session()
        check_active_runs(session, self.rpc_index)
        session.close()
