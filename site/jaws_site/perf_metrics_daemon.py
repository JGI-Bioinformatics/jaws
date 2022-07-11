import schedule
import time
import logging
from jaws_rpc import rpc_client_basic
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
        self.rpc_client = rpc_client_basic.RpcClientBasic(rpc_config, logger)

        # Set message expiration to 60 secs
        self.rpc_client.message_ttl = 3600

        self.perf_running_dir = config.conf.get("PERFORMANCE_METRICS", "running_dir")
        self.perf_done_dir = config.conf.get("PERFORMANCE_METRICS", "done_dir")
        self.perf_proc_dir = config.conf.get("PERFORMANCE_METRICS", "processed_dir")
        self.perf_cleanup_time = int(config.conf.get("PERFORMANCE_METRICS", "cleanup_time"))

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
        session = database.session_factory()
        performance_metrics = perf_metrics.PerformanceMetrics(session, self.rpc_client)
        performance_metrics.process_metrics(self.perf_done_dir, self.proc_dir)
        session.close()

    def cleanup(self):
        """
        Check for files in the running folder to move to done folder
        """

        # Get paths for both running/done
        running_dir = Path(self.perf_running_dir)
        done_dir = Path(self.perf_done_dir)

        # If there is no running folder return (there's probably an error here)
        if not running_dir.exists():
            logger.warning(f"Running folder not found: {running_dir}")
            return

        # If there is no done folder make sure to create one
        if not done_dir.exists():
            try:
                done_dir.mkdir(exist_ok=True)
            except Exception as ex:
                logger.warning(
                    f"Error making new directory {done_dir} {type(ex).__name__} : {ex}"
                )
        # Get all the csv files in the running dir
        files_running = running_dir.glob("*.csv")
        # Get the time the daemon was run
        now = time.time()
        for metric in files_running:
            try:
                status = metric.stat()
            except FileNotFoundError as ex:
                logger.warning(f"Error file not found {metric} : {ex}")
            except Exception as ex:
                logger.warning(f"Error getting file stat {type(ex).__name__} : {ex}")
            # gets the time the file hasn't been modified to in minutes
            idle_time = (now - status.st_ctime) / 60
            # If we're over the number of minutes and there has been no modifications then move it
            if idle_time > self.perf_cleanup_time:
                # Create a new path based on the metrics name and the done directory
                new_path = self.perf_done_dir / metric.name
                logger.debug(f"Moving {metric} to {new_path}")
                # Moves to the new path
                try:
                    metric.replace(new_path)
                except PermissionError as ex:
                    logger.warning(
                        f"Error moving file {metric} due to directory permissions {type(ex).__name__} : {ex}"
                    )
                except OSError as ex:
                    logger.warning(
                        f"Error moving file {metric} from one disk to another {type(ex).__name__} : {ex}"
                    )
