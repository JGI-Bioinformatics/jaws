"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import schedule
import time
import logging
from jaws_site import database, runs, config

# from jaws_site import perf_metrics_es
from jaws_rpc import rpc_client, rpc_client_basic

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
        self.central_rpc_client = rpc_client.RpcClient(
            config.conf.get_section("CENTRAL_RPC_CLIENT"), logger
        )
        self.report_rpc_client = rpc_client_basic.RpcClientBasic(
            config.conf.get_section("RUNS_ES_RPC_CLIENT"), logger
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
        session = database.Session()
        runs.check_active_runs(
            session, self.central_rpc_client, self.report_rpc_client
        )
        runs.send_run_status_logs(session, self.central_rpc_client)
        # perf_metrics_es.Metrics(session, self.pmetrics_es_rpc_client).process_metrics()
        session.close()
