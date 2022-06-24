import logging
import re
from typing import Callable
from pathlib import Path
from functools import lru_cache
from jaws_site import runs, tasks
from jaws_rpc import rpc_client_basic

import numpy as np
import pandas as pd

logger = logging.getLogger(__package__)


class PerformanceMetrics:
    """Run Performance Metrics"""
    def __init__(self, session: Callable, rpc_client: rpc_client_basic) -> None:
        self.session = session
        self.rpc_client = rpc_client
        self.tasks = {}

    @lru_cache()
    def get_run_id(self, cromwell_run_id: str) -> int:
        """Given a cromwell run id, return the jaws run id."""
        try:
            run = runs.Run.from_cromwell_run_id(self.session, cromwell_run_id)
        except runs.RunNotFoundError or runs.RunDbError as err:
            msg = f"Failed to get run_id from {cromwell_run_id=}: {err}"
            logger.warn(msg)
            raise runs.RunDbError(msg)
        return run.data.id

    def process_csv(self, csv_file: str) -> list:
        """Parse the performance csv files and generate a dictionary containing a list of performance metrics
        defined within the file.

        :param csv_file: performance metric csv file
        :type csv_file: str
        :return: list of dictionaries where each dictionary is a json doc of the performance metrics.
        :rtype: list
        """
        csv_data = pd.read_csv(csv_file, parse_dates=[0], index_col=[0])

        # Remove extranious parts from the current directory
        csv_data["current_dir"] = csv_data.current_dir.apply(remove_beginning_path)
        # Get data and make new columns in dataframe
        csv_data["cromwell_run_id"] = csv_data.current_dir.apply(extract_jaws_info)

        # Drops any non-workflow related processes
        csv_data = csv_data[csv_data.cromwell_run_id != "naw"]

        # If we have a cromwell-id then get the jaws id by calling jaws db
        csv_data["jaws_run_id"] = csv_data.cromwell_run_id.apply(self.get_run_id)
        # Filter out anything that didn't get a jaws_run_id returned
        csv_data = csv_data[csv_data.jaws_run_id != "None"]

        # Drops the current_dir, since we shouldn't need it anymore
        # csv_data = csv_data.drop(columns=["current_dir"])

        csv_data["@timestamp"] = csv_data.index.map(lambda x: x.isoformat())
        csv_data["mem_total"] = csv_data["mem_rss"] + csv_data["mem_vms"]
        csv_data["num_fds"] = csv_data["num_fds"].replace(["None"], np.nan)
        csv_data.fillna(0, inplace=True)
        csv_data["num_fds"] = csv_data.num_fds.astype(int)

        return csv_data.to_dict("records")

    @lru_cache()
    def add_taskname_mapping(self, cromwell_run_id):
        """Create a cache to map the cromwell task directory to the task name keyed using the cromwell_run_id.
        The cached result looks something like like:

        self.tasks[some_cromwell_run_id] = {
            "cromwell-executions/a/aa/aaa": "task_a",
            "cromwell-executions/b/bb/bbb": "task_b",
        }
        """
        try:
            self.tasks[cromwell_run_id] = tasks.TaskLog(
                self.session, cromwell_run_id=cromwell_run_id
            ).get_task_cromwell_dir_mapping()
        except Exception:
            # if calling TaskLog throws an exception (cromwell fails to get data using cromwell_run_id), don't
            # update the tasks list.
            pass

    def process_metrics(self, done_dir: str, proc_dir: str) -> None:
        """Gather all csv files in done_dir, parse and generate json docs for each entry, add docs to
        elasticsearch via RMQ submission whereby logstash picks up the doc from RMQ and inserts into
        the elasticsearch database.

        The csv files to be processed are stored in the done_dir. Once the file has been processed, it is moved
        to the proc_dir.

        :param done_dir: directory containing the performance csv files to process
        :type done_dir: str
        :param proc_dir: directory containing the processed performance csv files
        :type proc_dir: str
        :return: None
        """

        # If directories are not defined in config file, skip.
        if not done_dir or not proc_dir:
            return

        done_dir_obj = Path(done_dir)
        proc_dir_obj = Path(proc_dir)

        # If directories do not exists, skip.
        if not done_dir_obj.is_dir() or not proc_dir_obj.is_dir():
            return

        for done_file in list(done_dir_obj.glob("*.csv")):
            docs = self.process_csv(done_file)
            for doc in docs:
                cromwell_run_id = doc.get("cromwell_run_id")
                run_id = doc.get("jaws_run_id")

                # Just really a final check at this point since it should already be in the document
                if not run_id or not cromwell_run_id:
                    logger.error(
                        f"Error with {run_id=} or {cromwell_run_id=}, Not uploading to performance metrics"
                    )
                    continue

                logger.info(
                    f"Run {run_id}: Publish performance metrics for cromwell_run_id={cromwell_run_id}"
                )

                # Add cromwell dir to task name mapping if not already exists
                self.add_taskname_mapping(cromwell_run_id)

                # Add task name to doc
                cromwell_dir = doc["current_dir"]
                if (
                    cromwell_run_id in self.tasks
                    and cromwell_dir in self.tasks[cromwell_run_id]
                ):
                    task_name = self.tasks[cromwell_run_id][cromwell_dir]
                else:
                    task_name = "unknown"
                doc["task_name"] = task_name

                # Submit doc to RMQ to be picked up by logstash and inserted into elasticsearch
                self.rpc_client.request(doc)

            # Move csv file to processed folder
            processed_file = proc_dir_obj / done_file.name
            Path(f"{done_file}").rename(f"{processed_file}")

        # Clear cached cromwell dir to taskname mapping
        self.tasks = {}


@lru_cache()
def extract_jaws_info(working_dir):
    """
    Extract Cromwell run id from a Cromwell working directory
    returns "naw" if there is no cromwell id found

    :param task: command string
    :return: Cromwell run id in UUID format
    """
    cromwell_run_id = "naw"
    try:
        regex = re.compile(r"cromwell-executions\/[^\/]+\/([^\/]+)", re.I)
        match = regex.search(working_dir)
        if match:
            cromwell_run_id = match.group(1)
    except Exception as e:
        logger.warn(f"Error when processing {working_dir=}, {type(e).__name__} : {e}")

    return cromwell_run_id


@lru_cache()
def remove_beginning_path(working_dir):
    """Get the index of the string "cromwell-executions"""
    dir_name = "None"
    try:
        id_ce_dir = working_dir.find("cromwell-executions")
        # -1 means "cromwell-executions" was not found
        if id_ce_dir != -1:
            # Extract from "cromwell-executions" to the end of the string
            dir_name = working_dir[id_ce_dir:]
    except Exception as e:
        logger.warn(f"Error when processing {working_dir=}, {type(e).__name__} : {e}")

    return dir_name
