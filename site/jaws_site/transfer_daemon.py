import logging
import schedule
import time
from jaws_site import database, transfers


logger = logging.getLogger(__package__)


class TransferDaemon:
    """
    The Transfer Daemon periodically checks for new queued transfer tasks and submits them to the
    appropriate data transfer object (e.g. aws_s3_transfer).
    """

    def __init__(self):
        logger.info("Initializing transfer daemon")
        with database.session_factory() as session:
            transfers.reset_queue(session)

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
        with database.session_factory() as session:
            transfers.check_active_transfers(session)
