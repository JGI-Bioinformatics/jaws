"""
JAWS Daemon process periodically queries Globus to update xfer status and submit new xfers.
"""

import schedule
import time
import logging
from jaws_central.xfer_queue import XferQueue


logger = logging.getLogger(__package__)


class Daemon:
    """
    Daemon that periodically checks on this site's active runs.
    It prompts active runs to query Cromwell or Globus, as appropriate.
    When Globus uploads are completed, the run is submitted to Cromwell.
    When Cromwell execution has completed, the output is downloaded using Globus.
    When the download is finished, the run is complete.
    """

    def __init__(self):
        logger.info("Initializing daemon")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.check_active_xfers)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_active_xfers(self):
        """
        Query Globus for transfer task status updates and submit new xfers.
        """
        session = database.Session()
        xfer = XferQueue(session)
        xfer.update()
        session.close()
