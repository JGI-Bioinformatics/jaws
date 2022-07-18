import schedule
import time
import logging
from jaws_central.database import session_factory
from jaws_central.runs import check_active_runs


logger = logging.getLogger(__package__)


class RunDaemon:
    """
    Daemon that periodically checks on active runs.
    """

    def __init__(self):
        logger.info("Initializing run daemon")

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
        session = session_factory()
        check_active_runs(session)
        session.close()
