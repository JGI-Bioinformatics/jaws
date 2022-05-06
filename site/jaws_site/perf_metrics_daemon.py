import schedule
import time
import logging
import jaws_rpc
from jaws_site import config, database, perf_metrics

logger = logging.getLogger(__package__)


class PerformanceMetricsDaemon:
    """
    The Report Daemon periodically checks for newly completed Runs, generates the
    final Run report, and publishes to the elasticsearch service.
    """

    def __init__(self):
        logger.info("Initializing report daemon")
        rpc_config = config.conf.get_section("RUNS_ES_RPC_CLIENT")
        self.rpc_client = jaws_rpc.rpc_client_basic(rpc_config, logger)

        # Set message expiration to 60 secs
        self.runs_es_rpc_client.message_ttl = 3600
        self.pmetrics_es_rpc_client.message_ttl = 3600

        self.perf_done_dir = config.conf.get("PERFORMANCE_METRICS", "done_dir")
        self.perf_proc_dir = config.conf.get("PERFORMANCE_METRICS", "processed_dir")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(10).seconds.do(self.check_queue)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_queue(self):
        """
        Check for newly completed Runs.
        """
        session = database.Session()
        performance_metrics = perf_metrics.PerformanceMetrics(session, self.rpc_client)
        performance_metrics.process_metrics(self.perf_done_dir, self.proc_dir)
        session.close()
