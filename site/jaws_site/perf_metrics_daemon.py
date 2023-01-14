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
        rpc_config = config.conf.get_section("RMQ")
        rpc_config["queue"] = "JAWS_PERF_METRICS_ES"
        self.rpc_client = rpc_client_basic.RpcClientBasic(rpc_config, logger)

        # Set message expiration to 60 secs
        self.rpc_client.message_ttl = 3600

        self.running_dir = config.conf.get("PERFORMANCE_METRICS", "running_dir")
        self.done_dir = config.conf.get("PERFORMANCE_METRICS", "done_dir")
        self.proc_dir = config.conf.get("PERFORMANCE_METRICS", "processed_dir")
        self.error_dir = config.conf.get("PERFORMANCE_METRICS", "error_dir")
        self.perf_cleanup_time = int(config.conf.get("PERFORMANCE_METRICS", "cleanup_time"))

        # Create csv dirs if not exists
        done_path = Path(self.done_dir)
        proc_path = Path(self.proc_dir)
        error_path = Path(self.error_dir)
        try:
            done_path.mkdir(parents=True, exist_ok=True)
        except Exception as err:
            logger.error(f"Failed to create dir {done_path.as_posix()}: {err}")
            raise
        try:
            proc_path.mkdir(parents=True, exist_ok=True)
        except Exception as err:
            logger.error(f"Failed to create dir {proc_path.as_posix()}: {err}")
            raise
        try:
            error_path.mkdir(parents=True, exist_ok=True)
        except Exception as err:
            logger.error(f"Failed to create dir {error_path.as_posix()}: {err}")
            raise

    def start_daemon(self):
        """
        Run scheduled task(s) periodically.
        """
        schedule.every(self.perf_cleanup_time).seconds.do(self.cleanup)
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
        performance_metrics.process_metrics(done_dir=self.done_dir, proc_dir=self.proc_dir,
                                            error_dir=self.error_dir)
        session.close()

    def cleanup(self):
        """
        Check for files in the running folder to move to done folder
        """

        # Get paths for both running/done
        running_path = Path(self.running_dir)
        done_path = Path(self.done_dir)

        # If there is no running folder return (there's probably an error here)
        if not running_path.exists():
            logger.warning(f"Running folder not found: {running_path.as_posix()}")
            return

        # If there is no done folder make sure to create one
        if not done_path.exists():
            try:
                done_path.mkdir(exist_ok=True)
            except Exception as ex:
                logger.warning(
                    f"Error making new directory {done_path.as_posix()} {type(ex).__name__} : {ex}"  # noqa
                )
                return
        # Get all the csv files in the running dir
        files_running = running_path.glob("*.csv")
        # Get the time the daemon was run
        now = time.time()
        for metric in files_running:
            logger.debug(f"Checking to move {metric}")
            try:
                status = metric.stat()
            except FileNotFoundError as ex:
                logger.warning(f"Error file not found {metric} : {ex}")
                continue
            except Exception as ex:
                logger.warning(f"Error getting file stat {type(ex).__name__} : {ex}")
                continue

            # gets the time the file hasn't been modified to in minutes
            idle_time = (now - status.st_ctime) / 60
            # If we're over the number of minutes and there has been no modifications then move it
            if idle_time > int(self.perf_cleanup_time):
                try:
                    # Create a new path based on the metrics name and the done directory
                    new_path = done_path / metric.name
                    logger.info(f"Moving {metric} to {new_path}")
                    # Moves to the new path
                    metric.replace(new_path)
                except PermissionError as ex:
                    logger.warning(
                        f"Error moving file {metric} due to directory permissions {type(ex).__name__} : {ex}"
                    )
                    continue
                except OSError as ex:
                    logger.warning(
                        f"Error moving file {metric} from one disk to another {type(ex).__name__} : {ex}"
                    )
                    continue
                except Exception as ex:
                    logger.warning(
                        f"Error moving file {type(ex).__name__} : {ex}"
                    )
                    continue
