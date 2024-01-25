"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import logging
import time

import schedule
# from jaws_site import perf_metrics_es
from jaws_rpc import rpc_client, rpc_client_basic

from jaws_site import config, database, runs
from jaws_common.jaws_monitor import heartbeat

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
        central_rpc_params = config.conf.get_section("RMQ")
        central_rpc_params["queue"] = "CENTRAL"
        self.central_rpc_client = rpc_client.RpcClient(central_rpc_params, logger)
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

        # schedule the heartbeat for prometheus. The first arguement represents
        # the site, and the second, represents what the service is.
        site = self.central_rpc_params["queue"].lower()
        schedule.every(30).minutes.do(lambda: heartbeat(site, "run_daemon"))

        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_active_runs(self):
        """
        Check for runs in particular states.
        """
        with database.Session() as session:
            runs.check_active_runs(
                session, self.central_rpc_client, self.report_rpc_client
            )
            runs.send_run_status_logs(session, self.central_rpc_client)
