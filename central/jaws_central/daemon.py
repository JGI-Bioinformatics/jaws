"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import schedule
import time
import logging
from jaws_central import database, runs, config
from jaws_rpc import rpc_client


logger = logging.getLogger(__package__)


class Daemon:
    """
    Daemon that periodically checks on active Runs.
    For example, checks to see if a Run's files finished uploading successfully,
    then send the Run to the Site.  Otherwise, could send to a different Site.
    """

    def __init__(self):
        logger.info("Initializing daemon")
        site_rpc_params = config.conf.get_section("CENTRAL_RPC_CLIENT")
        self.rpc_index = rpc_index.rpc_index = rpc_index.RpcIndex(site_rpc_params, logger)

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
        runs.check_active_runs(session, self.rpc_index)
        session.close()
