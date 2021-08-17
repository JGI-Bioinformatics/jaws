"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import schedule
import time
import logging
from jaws_site import database, runs, config
from jaws_rpc import rpc_client


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
        self.central_rpc_client = rpc_client.RpcClient(config.conf.get_section("CENTRAL_RPC_CLIENT"), logger)

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
        runs.check_active_runs(session)
        runs.send_run_status_logs(session, self.central_rpc_client)
        session.close()
