import pandas as pd
import numpy as np
import logging
from typing import Callable
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from pathlib import Path
from jaws_site import config, models, rpc_es, runs_es

logger = logging.getLogger(__package__)


class RunDbError(Exception):
    pass


class Metrics:
    def __init__(self, session: Callable, rpc_client: rpc_es.RPCRequest) -> None:
        self.session = session
        self.rpc_client = rpc_client

    def get_run_id(self, cromwell_id: str) -> int:
        try:
            run = self.session.query(models.Run).filter(models.Run.cromwell_run_id == cromwell_id).one()
        except (IntegrityError, SQLAlchemyError) as err:
            msg = f"Failed to get run_id from {cromwell_id=}: {err}"
            logger.warn(msg)
            raise RunDbError(msg)
        return run.id

    def process_metrics(self):
        done_dir = config.conf.get("PERFORMANCE_METRICS", "done_dir")
        proc_dir = config.conf.get("PERFORMANCE_METRICS", "processed_dir")

        # If directories are not defined in config file, skip.
        if not done_dir or not proc_dir:
            return

        done_dir_obj = Path(done_dir)
        proc_dir_obj = Path(proc_dir)

        # If directories do not exists, skip.
        if not done_dir_obj.is_dir() or not proc_dir_obj.is_dir():
            return

        for done_file in list(done_dir_obj.glob('*.csv')):
            docs = process_csv(done_file)
            for doc in docs:
                cromwell_id = doc.get('cromwell_id')
                if not cromwell_id:
                    logger.warn(f"Publishing performance metrics: cannot find cromwell_run_id in {done_file=}.")
                    continue
                try:
                    run_id = self.get_run_id(cromwell_id)
                except RunDbError as err:
                    logger.warn(f"Publishing performace metrics: failed to get run_id from {cromwell_id=}: {err})")
                    continue
                doc['jaws_run_id'] = run_id
                logger.info(f"Run {run_id}: Publish performance metrics for cromwell_id={cromwell_id}")
                response, status = runs_es.send_rpc_run_metadata(self.rpc_client, doc)

            # Move csv file to processed folder
            processed_file = proc_dir_obj / done_file.name
            Path(f"{done_file}").rename(f"{processed_file}")


def extract_jaws_info(working_dir):
    try:
        # Make defaults in case there are no values later
        default_response = "naw", "naw", "naw", "naw", \
            "naw", "naw", 0  # not a workflow
        # no subworflow
        sub_workflow_name = sub_cromwell_id = sub_task_name = "nosub"
        shard_n = 0

        # Make dummy varialbes
        workflow_name = None
        task_name = None
        cromwell_id = None

        # Splits the directory structure into a list
        try:
            split = working_dir.split("/")
        except AttributeError:
            return default_response

        # If there is no "cromwell-executions" in the directory structure
        # then it's not a part of cromwell executions
        # and therefore not a part of the JAWS workflow
        # Label it so that it can be removed later
        if "cromwell-executions" not in split:
            return default_response

        # Only keep the parts after "cromwell-executions"
        split = split[split.index("cromwell-executions")+1:]
        split.remove("execution")

        # Fill in outputs from the split string by looping through the parts
        for s in split:
            # Get task information
            if "call-" in s:
                # If we don't have a task yet it's the main task
                if task_name is None:
                    task_name = s.split('-')[-1]
                # Otherwise it's a subtask
                else:
                    sub_task_name = s.split('-')[-1]
            # Get information for shards
            elif "shard-" in s:
                shard_n = s.split('-')[-1]
            # Cromwell-id should be 8-4-4-4-12 giving 5 parts
            # If we can split it into 5 parts it's a {sub-}workflow
            elif len(s.split("-")) == 5:
                if cromwell_id is None:
                    cromwell_id = s
                else:
                    sub_cromwell_id = s
            # Any other strings should be the {sub-}workflow name
            else:
                if workflow_name is None:
                    workflow_name = s
                else:
                    sub_workflow_name = s

        return workflow_name, cromwell_id, task_name, sub_workflow_name, \
            sub_cromwell_id, sub_task_name, shard_n
    except Exception as e:
        logger.warn(f"Error when processing cromwell_id={cromwell_id}, {type(e).__name__} : {e}")
        return default_response


def process_csv(csv_file):
    csv_data = pd.read_csv(csv_file, parse_dates=[0], index_col=[0])

    # Get data and make new columns in dataframe
    csv_data["workflow_name"], csv_data["cromwell_id"], \
        csv_data["task_name"], csv_data["sub_workflow_name"], \
        csv_data["sub_cromwell_id"], csv_data["sub_task_name"], \
        csv_data["shard_n"] = csv_data.current_dir.apply(extract_jaws_info)
    # Drops any non-workflow related processes
    csv_data = csv_data[csv_data.workflow_name != "naw"]
    # Drops the current_dir, since we shouldn't need it anymore
    csv_data = csv_data.drop(columns=['current_dir'])

    csv_data['@timestamp'] = csv_data.index.map(lambda x: x.isoformat())
    csv_data['mem_total'] = csv_data['mem_rss'] + csv_data['mem_vms']
    csv_data['num_fds'] = csv_data['num_fds'].replace(['None'], np.nan)
    csv_data.fillna(0, inplace=True)
    csv_data["num_fds"] = csv_data.num_fds.astype(int)

    return csv_data.to_dict('records')
