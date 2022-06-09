import logging
import numpy as np
import datetime
from typing import Tuple, Any
from jaws_central import dateutils, es_operations
from jaws_central.es_operations import ElasticsearchError

logger = logging.getLogger(__package__)


class Report:
    def __init__(self, es_client: es_operations.ESClient) -> None:
        self.es_client = es_client

    @staticmethod
    def parse_perf_metrics(pmetrics: list, datetime_format: str) -> Tuple[list, list, str]:
        """Given a list of performance metric docs from elasticsearch in dict format, gather all cpus, memory,
        num_threads into lists. All docs pertain to a specific jaws task. Also return start and end date for the
        task.

        :param pmetrics: list of dicts containing performance metrics obtained from elasticsearch
        :type pmetrics: list
        :param datetime_format: date format of the jaws task
        :type datetime_format: str
        :return a tuple of lists containing cpus, memory, threads, start and end dates
        """
        cpus = []
        mems = []
        threads = []
        start_dt = None
        end_dt = None

        for pmetric in pmetrics:
            if not pmetric:
                continue
            cpu_user = pmetric.get("cpu_user", 0)
            cpu_system = pmetric.get("cpu_system", 0)
            num_threads = pmetric.get("num_threads", 0)
            mem = pmetric.get("mem_total", 0)
            cpus.append((cpu_user + cpu_system) / num_threads)
            mems.append(mem)
            threads.append(num_threads)

            timestamp = dateutils.strip_timezone_from_date(pmetric["@timestamp"])
            task_dt = datetime.datetime.strptime(timestamp, datetime_format)
            if not start_dt or task_dt < start_dt:
                start_dt = task_dt
            if not end_dt or task_dt > end_dt:
                end_dt = task_dt

        return cpus, mems, threads, start_dt, end_dt

    def get_task_perf_metrics(self, pmetrics: list, datetime_format: str) -> Tuple[dict, float]:
        """Given a list of performance metric docs from elasticsearch in dict format, compute min, mean, max
        for cpu, memory, num threads and runtime. The performance metric docs correspond to time series data
        points for a given jaws task.

        :param pmetrics: list of dicts containing performance metrics obtained from elasticsearch
        :type pmetrics: list
        :param datetime_format: date format of the jaws task
        :type datetime_format: str
        :return a tuple of dict containing the performance metric key/value pairs and the runtime in secs.
        """

        cpus, mems, threads, start_dt, end_dt = self.parse_perf_metrics(pmetrics, datetime_format)
        start_time = datetime.datetime.strftime(start_dt, datetime_format)
        end_time = datetime.datetime.strftime(end_dt, datetime_format)
        runtime_sec = dateutils.get_seconds_between_dates(start_time, end_time, datetime_format)

        cpu_array = np.array(cpus)
        mem_array = np.array(mems)

        decimal_precision = 2
        cpu_mean = float(round(np.mean(cpu_array), decimal_precision))
        cpu_min = float(round(np.min(cpu_array), decimal_precision))
        cpu_max = float(round(np.max(cpu_array), decimal_precision))
        mem_mean = float(round(np.mean(mem_array), decimal_precision))
        mem_min = float(round(np.min(mem_array), decimal_precision))
        mem_max = float(round(np.max(mem_array), decimal_precision))
        threads_mean = float(round(np.mean(threads), decimal_precision))
        threads_min = float(round(np.min(threads), decimal_precision))
        threads_max = float(round(np.max(threads), decimal_precision))
        cpu_mean_pct = float(round(cpu_mean / runtime_sec * 100, decimal_precision)) if runtime_sec else 0
        cpu_min_pct = float(round(cpu_min / runtime_sec * 100, decimal_precision)) if runtime_sec else 0
        cpu_max_pct = float(round(cpu_max / runtime_sec * 100, decimal_precision)) if runtime_sec else 0

        metrics = {}
        metrics["cpu_mean"] = cpu_mean
        metrics["cpu_min"] = cpu_min
        metrics["cpu_max"] = cpu_max
        metrics["cpu_mean_pct"] = cpu_mean_pct
        metrics["cpu_min_pct"] = cpu_min_pct
        metrics["cpu_max_pct"] = cpu_max_pct
        metrics["memory_mean"] = mem_mean
        metrics["memory_min"] = mem_min
        metrics["memory_max"] = mem_max
        metrics["num_threads_mean"] = threads_mean
        metrics["num_threads_min"] = threads_min
        metrics["num_threads_max"] = threads_max
        metrics["runtime_sec"] = runtime_sec

        return metrics, runtime_sec

    def add_perf_metrics_to_run_metadata(self, run_metadata: dict) -> dict:
        """Given a run metadata doc (dict) from elasticsearch, lookup the performance metrics for given run and
        add them to the dict. Return dict with added performance metrics.

        :param run_metadata: dict containing run metadata from elasticsearch
        :type run_metadata: dict
        :return dict
        """
        run_id = run_metadata["run_id"]
        tasks = run_metadata["tasks"]
        task_datetime_format = "%Y-%m-%d %H:%M:%S"
        pmetric_datetime_format = "%Y-%m-%dT%H:%M:%S"
        submit_time = run_metadata["submitted"]

        new_tasks = []
        total_runtime_sec = 0
        total_max_cpu = 0
        total_max_mem = 0
        total_max_threads = 0
        run_max_cpu_pct = []
        run_min_cpu_pct = []
        run_mean_cpu_pct = []
        pmetrics_found = False

        for task in tasks:
            task_name = task["name"]

            # Each task may contain multiple pmetric datapoints. Here, we're gathering all performance metrics
            # for this task.
            pmetrics = self.es_client.search_perf_metrics_by_task(run_id, task_name)

            if not pmetrics:
                continue

            # If performance metrics found for a task, set flag
            pmetrics_found = True

            # Get statistics of compute metrics for this run task
            task_pmetrics, runtime_sec = self.get_task_perf_metrics(pmetrics, pmetric_datetime_format)

            if "run_start" in task:
                task_start_time = task["run_start"]
                wait_time_sec = dateutils.get_seconds_between_dates(submit_time, task_start_time, task_datetime_format)
            else:
                # Set wait time to -1 if this value can't be determined
                wait_time_sec = -1

            task_pmetrics["wait_time_sec"] = wait_time_sec

            # Add performance metric stats to run task metadada
            task.update(task_pmetrics)

            new_tasks.append(task)

            total_runtime_sec += runtime_sec

            max_mem = task.get("memory_max", 0)
            if max_mem > total_max_mem:
                total_max_mem = max_mem

            max_threads = task.get("num_threads_max", 0)
            if max_threads > total_max_threads:
                total_max_threads = max_threads

            total_max_cpu += task["cpu_max"]

            run_min_cpu_pct.append(task["cpu_min_pct"])
            run_max_cpu_pct.append(task["cpu_max_pct"])
            run_mean_cpu_pct.append(task["cpu_mean_pct"])

        # If no performance metrics found for any task, return None to indicate that no metrics were found.
        if not pmetrics_found:
            return None

        # del run_metadata["@version"]

        total_start_date = run_metadata["submitted"]
        total_end_date = run_metadata["updated"]

        run_metadata["total_walltime_sec"] = dateutils.get_seconds_between_dates(total_start_date, total_end_date,
                                                                                 "%Y-%m-%d %H:%M:%S")
        run_metadata["total_runtime_sec"] = total_runtime_sec
        run_metadata["cpu_max"] = total_max_cpu
        run_metadata["memory_max"] = total_max_mem
        run_metadata["num_threads_max"] = total_max_threads
        run_metadata["tasks"] = new_tasks

        run_min_cpu_pct_array = np.array(run_min_cpu_pct)
        run_max_cpu_pct_array = np.array(run_max_cpu_pct)
        run_mean_cpu_pct_array = np.array(run_mean_cpu_pct)

        run_metadata["cpu_max_pct"] = np.mean(run_max_cpu_pct_array)
        run_metadata["cpu_min_pct"] = np.mean(run_min_cpu_pct_array)
        run_metadata["cpu_mean_pct"] = np.mean(run_mean_cpu_pct_array)

        return run_metadata

    def update_run_metadata_with_perf_metrics(self, run_id: int) -> Any:
        """Given a jaws run id, lookup run metadata doc from elasticsearch, merge performance metrics to doc
        and update doc in elasticsearch.

        :param run_id: jaws run id to update
        :type run_id: int
        :return elasticsearch response from inserting update doc into elasticsearch
        """
        entry = self.es_client.search_by_run_id(run_id)
        if not entry:
            return

        run_metadata = entry[0]

        try:
            new_entry = self.add_perf_metrics_to_run_metadata(run_metadata)
        except ElasticsearchError:
            raise

        if new_entry:
            try:
                response = self.es_client.insert_into_elasticsearch(new_entry)
            except ElasticsearchError:
                raise
        else:
            response = {"no_performance_entries_found": True}

        return response

    def get_run_ids_without_perf_metrics(self, start_date, end_date) -> list:
        """Given a start and end date string in the format YYYY-MM-DD, lookup all jaws run metadata in
        elasticsearch that doesn't have performance metrics. Return a list of run ids without those metrics.

        :param start_date: starting date to query elasticsearch for run metadata
        :type start_date: str in format YYYY-MM-DD
        :param end_date: ending date to query elasticsearch for run metadata
        :type end_date: str in format YYYY-MM-DD
        :return list of run ids
        """
        try:
            docs = self.es_client.search_run_by_dates(start_date=start_date, end_date=end_date,
                                                      has_no_perf_metrics_only=True)
        except ElasticsearchError:
            raise

        run_ids = []
        for doc in docs:
            if "run_id" in doc:
                run_ids.append(doc["run_id"])

        return run_ids
