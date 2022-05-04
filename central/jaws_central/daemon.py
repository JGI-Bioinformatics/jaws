import schedule
import time
import logging
from jaws_central import database
from jaws_central.runs import check_active_runs


logger = logging.getLogger(__package__)


class Daemon:
    """
    Daemon that periodically checks on this site's active runs.
    It prompts active runs to query Cromwell or Globus, as appropriate.
    When Globus uploads are completed, the run is submitted to Cromwell.
    When Cromwell execution has completed, the output is downloaded using Globus.
    When the download is finished, the run is complete.
    """

    def __init__(self, rpc_index):
        logger.info("Initializing daemon")
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
