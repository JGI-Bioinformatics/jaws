import schedule
import time
import logging
from jaws_central import dateutils, run_report, es_operations


logger = logging.getLogger(__package__)


class ReportDaemon:
    """Daemon that periodically updates run metadata in elasticsearch with performance metrics."""

    def __init__(self) -> None:
        logger.info("Initializing run daemon")
        try:
            self.es_client = es_operations.ESClient()
        except es_operations.ElasticsearchError as err:
            logger.error(err)
            raise
        self.report_num_days = 7
        self.schedule_time_sec = 60

    def start_daemon(self) -> None:
        """Run scheduled task(s) periodically."""
        schedule.every(self.schedule_time_sec).seconds.do(self.update_run_metadata)
        while True:
            schedule.run_pending()
            time.sleep(1)

    def update_run_metadata(self) -> None:
        """Lookup run metadata in elasticsearch that doesn't have performance metrics and
        update metadata with corresponding metrics based on a date range.
        """
        curr_date = dateutils.today()
        start_date = dateutils.subtract_days(curr_date, num_days=self.report_num_days)
        run_ids = run_report.Report(self.es_client).get_run_ids_without_perf_metrics(start_date, curr_date)

        for run_id in run_ids:
            report = run_report.Report(self.es_client)
            response = report.update_run_metadata_with_perf_metrics(run_id)
            logger.debug(f"Updating {run_id=}: {response=}")
