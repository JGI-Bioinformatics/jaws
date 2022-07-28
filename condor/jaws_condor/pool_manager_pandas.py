import logging

from typing import Dict

import pandas as pd
from jaws_condor.cmd_utils import run_sh_command
from jaws_condor.htcondor_cmds import HTCondor
import math

logger = logging.getLogger(__package__)

# TODO : These should all be created from a configuration file at some point
user_name = "svc-jtm-user"
condor_root = ".."
compute_types = ["medium", "xlarge"]

wanted_columns = "ClusterId RequestMemory RequestCpus \
    CumulativeRemoteSysCpu CumulativeRemoteUserCpu JobStatus \
        NumShadowStarts JobRunCount RemoteHost JobStartDate QDate"
condor_q_cmd = f"condor_q -allusers -af {wanted_columns}"
condor_idle_nodes = 'condor_status -const "TotalSlots == 1" -af Machine'

COMPUTE_SITE = "tahoma"
MIN_POOL = {}
MIN_POOL['medium'] = 4
MIN_POOL['xlarge'] = 0
MAX_POOL = 30

worker_sizes_sites = {
    "nersc": {
        "medium_cpu": 64,
        "medium_mem": 120,
        "xlarge_cpu": 72,
        "xlarge_mem": 1500,
    },
    "tahoma": {
        "medium_cpu": 72,
        "medium_mem": 364,
        "xlarge_cpu": 72,
        "xlarge_mem": 1480,
    }
}

cpu_bins = {}
cpu_labels = {}
cpu_bins['nersc'] = [0, 64, 72, 10_000_000]
cpu_labels['nersc'] = ["medium", "xlarge", "over"]

mem_bins = {}
mem_labels = {}
mem_bins['nersc'] = [0, 120, 1500, 10_000_000]
mem_labels['nersc'] = ["medium", "xlarge", "over"]

worker_sizes = worker_sizes_sites[COMPUTE_SITE]

squeue_args = {}
squeue_args["nersc"] = "--noheader --clusters=all -p genepool,genepool_shared,exvivo,exvivo_shared"
squeue_args["tahoma"] = "--noheader"
squeue_cmd = f'squeue -u {user_name} --format="%.12i %.10P %.50j %.10T %.10L %.12R" {squeue_args[COMPUTE_SITE]}'
squeue_columns = ["JOBID", "PARTITION", "NAME", "STATE", "TIME_LEFT", "NODELIST"]
# TODO : These should all be created from a configuration file at some point


def slurm_time_to_sec(time_str):
    # split off days first
    time_str = time_str.split("-")
    if len(time_str) > 1:
        days = int(time_str[0])
    else:
        days = 0

    # split time_str into HH:MM:SS
    time_str = time_str[-1].split(":")

    time_str_bits = {0: 1, 1: 60, 2: 60*60, 3: 60*60*24}
    total = 0
    # Run in reverse becasue we will always
    # have seconds and not always hours
    # 0 -> sec
    # 1 -> min
    # 2 -> hrs.
    # 3 -> days.
    for i, t in enumerate(time_str[::-1]):
        total += (time_str_bits[i]*int(t))

    if days > 0:
        total += (time_str_bits[3]*int(days))

    # Return total seconds
    return total


class PoolManagerPandas:
    """Class representing a single Run"""

    def __init__(self, condor_provider=None, slurm_provider=None, **kwargs):
        self.site = "nersc"
        self.condor_provider = condor_provider
        self.slurm_provider = slurm_provider
        # config.conf.get("SITE", "id"),

    def add_workers():
        logger.info("Checking to add workers to pool")
        return None

    def rm_workers():
        logger.info("Checking to remove workers to pool")
        return None

    def get_current_slurm_workers(self) -> Dict:
        slurm_status = {}

        _stdout, _stderr, _errcode = run_sh_command(squeue_cmd, show_stdout=False)

        # NOTE: Leaving in for furute SuperFacility calls at nersc
        # _stdout = sfapi.custom_cmd(
        #     token=access_token.token, cmd=squeue_cmd, site=site)
        # try:
        #     _stdout = _stdout['output']
        # except KeyError:
        #     logger.error("No output from squeue")
        #     output = {}
        #     for _type in compute_types:
        #         for _state in ["running", "pending"]:
        #             output[f'{_type}_{_state}']
        #     return output

        # If we have an error return a dictionary with 0 for each type and state
        if _errcode != 0:
            logger.error(f"No output from squeue: {squeue_cmd}")
            output = {}
            for _type in compute_types:
                for _state in ["running", "pending"]:
                    output[f'{_type}_{_state}'] = 0
            return output

        # Gets jobs from output by splitting on new lines
        jobs = [job.split() for job in _stdout.split("\n")]

        # Make dataframe to query with
        df = pd.DataFrame(jobs, columns=squeue_columns)
        # Drops rows if they have nan values
        df = df.dropna(axis=0)
        df['TIME_SEC'] = df["TIME_LEFT"].apply(slurm_time_to_sec)

        # Selections for running and pending
        mask_pending = (df["STATE"] == "PENDING")
        mask_running = (df["STATE"] == "RUNNING")

        # Selections for just condor jobs
        mask_condor = df["NAME"].str.contains("condor")

        # For sites with multiple types we want to know
        for _type in compute_types:
            mask_type = df["NAME"].str.contains(_type)
            # Each of these selects for a certian type of node based on a set of masks
            # Add the number of nodes to get how many are in each catogory
            slurm_status[f"{_type}_pending"] = sum(mask_type & mask_condor & mask_pending)
            slurm_status[f"{_type}_running"] = sum(mask_type & mask_condor & mask_running)

        slurm_running_df = df[mask_condor]
        logger.info(slurm_status)
        return slurm_status, slurm_running_df

    def get_condor_job_queue(self) -> pd.DataFrame:
        condor_jobs = self.condor_provider.condor_q()
        return condor_jobs

    def determine_condor_job_sizes(self, condor_jobs):
        df = pd.DataFrame(condor_jobs)
        df["mem_bin"] = pd.cut(
            df["RequestMemory"],
            bins=mem_bins[COMPUTE_SITE],
            labels=mem_labels[COMPUTE_SITE],
        )

        df["cpu_bin"] = pd.cut(
            df["RequestCpus"],
            bins=cpu_bins[COMPUTE_SITE],
            labels=cpu_labels[COMPUTE_SITE],
        )

        condor_q_status = {}
        mask_running_status = df["JobStatus"].astype(int) == 2
        mask_idle_status = df["JobStatus"].astype(int) == 1
        mask_hold_status = df["JobStatus"].astype(int) == 5

        mask_over = df["cpu_bin"].str.contains(
            "over") | df["mem_bin"].str.contains("over")

        condor_q_status["hold_and_impossible"] = sum(mask_over | mask_hold_status)

        for _type in compute_types:
            mask_mem_type = df["mem_bin"].str.contains(f"{_type}")
            # mask_cpu_type = df["cpu_bin"].str.contains(f"{_type}")
            mask_type = (mask_mem_type & ~(mask_over | mask_hold_status))

            condor_q_status[f"idle_{_type}"] = sum(
                mask_type & mask_idle_status
            )
            condor_q_status[f"running_{_type}"] = sum(
                mask_type & mask_running_status
            )

            condor_q_status[f"{_type}_cpu_needed"] = sum(
                df[mask_type].RequestCpus)
            condor_q_status[f"{_type}_mem_needed"] = sum(
                df[mask_type].RequestMemory)

        return condor_q_status

    def need_new_nodes(self, condor_job_queue: Dict, slurm_workers: Dict, machine_size: str) -> Dict:
        """
        Using the two dictionaries from the condor_q and squeue
        determine if we need any new workers for the machine_size types.
        """
        workers_needed = 0

        # Determines how many full (or partially full nodes) we need to create
        _cpu = (
            condor_job_queue[f"{machine_size}_cpu_needed"] /
            worker_sizes[f"{machine_size}_cpu"]
        )
        _mem = (
            condor_job_queue[f"{machine_size}_mem_needed"] /
            worker_sizes[f"{machine_size}_mem"]
        )
        _cpu = math.ceil(_cpu)
        _mem = math.ceil(_mem)

        workers_needed += max(_cpu, _mem)

        # If full workers_needed is 0 but we have work to be done still get a node
        if workers_needed == 0:
            if condor_job_queue[f"{machine_size}_cpu_needed"] or condor_job_queue[f"{machine_size}_mem_needed"]:
                workers_needed = 1

        # Total number running and pending to run (i.e. worker pool)
        current_pool_size = (
            slurm_workers[f"{machine_size}_pending"]
            + slurm_workers[f"{machine_size}_running"]
        )

        # If workers_needed is higher than the pool we'll add the diference
        # Else we don't need workers (add 0)
        # workers_needed = max(0, workers_needed - current_pool_size)
        workers_needed = (workers_needed - current_pool_size)
        if workers_needed < 0:
            workers_needed = (abs(workers_needed) - current_pool_size)  # + abs(workers_needed)

        if abs(workers_needed) < MIN_POOL[machine_size]:
            return 0
        # If we have less running than the minimum we always need to add more
        # Either add what we need from queue (workers_needed)
        # Or what we're lacking in the pool (min - worker pool)

        if current_pool_size < MIN_POOL[machine_size]:
            workers_needed = max(MIN_POOL[machine_size] - current_pool_size, workers_needed)

        # Check to make sure we don't go over the max pool size
        if (workers_needed + current_pool_size) > MAX_POOL:
            # Only add up to max pool and no more
            workers_needed = MAX_POOL - current_pool_size

        return workers_needed if workers_needed > 0 else 0

    def need_cleanup(self, condor_job_queue: Dict, slurm_workers: Dict, machine_size: str) -> Dict:
        """
        Using the two dictionaries from the condor_q and squeue
        determine if we need any new workers for the machine_size types.
        """
        # Start off with minimum pool
        workers_needed = MIN_POOL[machine_size]

        # Determines how many full (or partially full nodes) we need to create
        _cpu = (
            condor_job_queue[f"{machine_size}_cpu_needed"] /
            worker_sizes[f"{machine_size}_cpu"]
        )
        _mem = (
            condor_job_queue[f"{machine_size}_mem_needed"] /
            worker_sizes[f"{machine_size}_mem"]
        )

        # Round the numbers up
        _cpu = math.floor(_cpu)
        _mem = math.floor(_mem)

        workers_needed += max(_cpu, _mem)

        # Total number running and pending to run (i.e. worker pool)
        current_pool_size = (
            slurm_workers[f"{machine_size}_pending"]
            + slurm_workers[f"{machine_size}_running"]
        )

        # If workers_needed is higher than the pool we'll add the diference
        # Else we don't need workers (add 0)
        workers_needed = (workers_needed - current_pool_size)

        if workers_needed >= MIN_POOL[machine_size]:
            workers_needed = workers_needed - MIN_POOL[machine_size]

        return workers_needed if workers_needed < 0 else 0

    def run_cleanup(self, slurm_running_df, cleanup_num, _type):
        # Runs a condor_q autoformat to get the desired columns back

        # Gets the idle nodes from condor
        idle_nodes = self.condor_provider.condor_idle()
        logger.info(idle_nodes)

        try:
            logger.info(slurm_running_df.shape)
        except NameError as e:
            logger.info(f'No slurm nodes yet, {e}')
            return None

        slurm_running_df.sort_values("TIME_SEC", inplace=True, ascending=False)

        # Any pending nodes are at the front of the list to remove
        pending = slurm_running_df[slurm_running_df.STATE == "PENDING"].NODELIST
        nodes = [n for n in pending]

        # Add nodes that are idle to the end of the list
        nodes.extend(idle_nodes)

        # If there are no nodes to remove just return
        if len(nodes) == 0:
            return None

        num = 0
        for node in nodes:
            if num >= cleanup_num:
                continue
            try:
                node_mask = (slurm_running_df['NODELIST'] == node)
                type_mask = slurm_running_df["NAME"].str.contains(_type)
                job_id = slurm_running_df[node_mask & type_mask]
                job_id = int(job_id.JOBID)
            except IndexError:
                continue
            except TypeError:
                continue
            num += 1
            logger.info(f"Removing {node} with JobID {job_id}")
            # print(f"{job_id} ", end="")

            # _stdout, _stderr, _errorcode = run_sh_command(f"scancel {job_id}", show_stdout=False)
            # if _errorcode != 0:
            #     print(f"ERROR: failed to execute condor_q command: scancel {job_id} : {_stderr}")
            #     return None

            # x = sfapi.delete_job(access_token.token,
            #                      site=COMPUTE_SITE, jobid=job_id)
            # logger.info(x)

        return nodes

    def run_sbatch(self, new_workers: int = 0, _type: str = "medium"):
        for i in range(new_workers):
            print(f"Running sbatch {_type}")


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    pool = PoolManagerPandas(condor_provider=HTCondor())
    slurm_status, slurm_running_df = pool.get_current_slurm_workers()
    print(slurm_status)
    print(slurm_running_df)
    condor_status = pool.get_condor_job_queue()
    print(condor_status)
    work_status = pool.determine_condor_job_sizes(condor_status)
    print(work_status)
    for _type in compute_types:
        old_workers = pool.need_cleanup(work_status, slurm_status, _type)
        new_workers = pool.need_new_nodes(work_status, slurm_status, _type)
        print(f"{old_workers}\t{new_workers}")
        if old_workers < 0:
            pool.run_cleanup(slurm_running_df, abs(old_workers), _type)
        if new_workers > 0:
            pool.run_sbatch(abs(new_workers), _type)
