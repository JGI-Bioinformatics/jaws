import schedule
import time
import logging
import jaws_rpc
from jaws_site import config, database, perf_metrics
from pathlib import Path

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

        self.perf_running_dir = config.conf.get("PERFORMANCE_METRICS", "running_dir")
        self.perf_done_dir = config.conf.get("PERFORMANCE_METRICS", "done_dir")
        self.perf_proc_dir = config.conf.get("PERFORMANCE_METRICS", "processed_dir")
        self.perf_cleanup_time = config.conf.get("PERFORMANCE_METRICS", "cleanup_time")

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(self.perf_cleanup_time).minutes.do(self.cleanup)
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

    def cleanup(self):
        """
        Check for files in the running folder to move to done folder
        """
        try:
            # Get paths for both running/done
            running_dir = Path(self.perf_running_dir)
            done_dir = Path(self.perf_done_dir)

            # If there is no running folder return (there's probably an error here)
            if not running_dir.exists():
                logger.warn(f"Running folder not found: {running_dir}")
                return

            # If there is no done folder make sure to create one
            if not done_dir.exists():
                done_dir.mkdir(exist_ok=True)

            # Get all the csv files in the running dir
            files_running = running_dir.glob("*.csv")
            # Get the time the daemon was run
            now = time.time()
            for metric in files_running:
                status = metric.stat()
                # gets the time the file hasn't been modified to in minutes
                idle_time = (now-status.st_mtime)/60
                # If we're over the number of minutes and there has been no modifications then move it
                if idle_time > self.perf_cleanup_time:
                    # Create a new path based on the metrics name and the done directory
                    new_path = self.perf_done_dir / metric.name
                    logger.debug(f"Moving {metric} to {new_path}")
                    # Moves to the new path
                    metric.replace(new_path)
        except Exception as ex:
            logger.warn(f"Unknown exception {type(ex).__name__} : {ex}")
