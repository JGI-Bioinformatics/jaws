import logging
import time

import schedule

from jaws_site import transfers

logger = logging.getLogger(__package__)


def write_heartbeat():
    logger.info("Transfer daemon heartbeat: it's up and running")


class TransferDaemon:
    """
    The Transfer Daemon periodically checks for new queued transfer tasks and submits them to the
    appropriate data transfer object (e.g. aws_s3_transfer).
    """

    def __init__(self):
        logger.info("Initializing transfer daemon")
        # clean up any interrupted transfers so they can be retried
        transfers.reset_queue()

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(30).minutes.do(write_heartbeat)
        schedule.every(10).seconds.do(transfers.check_fix_perms_queue)
        schedule.every(10).seconds.do(transfers.check_transfer_queue)
        while True:
            schedule.run_pending()
            time.sleep(1)
