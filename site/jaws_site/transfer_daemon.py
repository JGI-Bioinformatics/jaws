import schedule
import time
import logging
from jaws_site import database, transfers


logger = logging.getLogger(__package__)


class TransferDaemon:
    """
    The Transfer Daemon periodically checks for new queued transfer tasks and submits them to the
    appropriate data transfer object (e.g. aws_s3_transfer).
    """

    def __init__(self):
        logger.info("Initializing transfer daemon")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.check_transfer_queue)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_transfer_queue(self):
        """
        Do any queued transfers now.
        """
        session = database.Session()
        transfers.check_queue(session)
        session.close()