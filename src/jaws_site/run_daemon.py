"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import logging
import time

import schedule

# from jaws_site import perf_metrics_es
from jaws_rpc import rpc_client_basic

from jaws_site import config, database, runs

logger = logging.getLogger(__package__)


class RunDaemon:
    """
    Daemon that periodically checks on this site's active runs.
    It prompts active runs to query Cromwell or Globus, as appropriate.
    When Globus uploads are completed, the run is submitted to Cromwell.
    When Cromwell execution has completed, the output is downloaded using Globus.
    When the download is finished, the run is complete.
    """

    def __init__(self):
        logger.info("Initializing daemon")
        report_rpc_params = config.conf.get_section("RMQ")
        report_rpc_params["queue"] = "RUNS_ES_RPC_CLIENT"
        self.report_rpc_client = rpc_client_basic.RpcClientBasic(
            report_rpc_params, logger
        )

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
        with database.Session() as session:
            runs.check_active_runs(session, self.report_rpc_client)
            runs.send_run_logs(session)
