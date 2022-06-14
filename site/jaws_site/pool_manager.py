import logging
from typing import Dict

from jaws_site import config
import pandas as pd
import time
from jaws_site.utils import run_sh_command

logger = logging.getLogger(__package__)

# System vars, dirs, and cmds
#
user_name = "jaws_jtm"
condor_root = ".."
compute_types = ["regular", "xlarge"]

wanted_columns = "ClusterId RequestMemory RequestCpus CumulativeRemoteSysCpu CumulativeRemoteUserCpu JobStatus NumJobStarts RemoteHost JobStartDate QDate"
condor_q_cmd = f"condor_q -allusers -af {wanted_columns}"
condor_idle_nodes = 'condor_status -const "TotalSlots == 1" -af Machine'


###################### TODO : These should all be created from a configuration file at some point ######################
COMPUTE_SITE = "cori"
MIN_POOL = 0
MAX_POOL = 10
MAX_SUBMIT_SIZE = 10

worker_sizes_sites = {"perlmutter": {
    "regular_cpu": 256,
    "regular_mem": 500,
},
    "cori": {
    "regular_cpu": 64,
    "regular_mem": 120,
}}

worker_sizes = worker_sizes_sites[COMPUTE_SITE]

# Extra partition arguments for cori
squeue_args = {}
squeue_args["cori"] = "--clusters=all -p genepool,genepool_shared,exvivo,exvivo_shared"
squeue_args["perlmutter"] = "--clusters=all"
squeue_cmd = f'squeue -u {user_name} --format="%.18i %.24P %.100j %.10T %S %e %.70R" {squeue_args[COMPUTE_SITE]}'
squeue_columns = ["JOBID", "PARTITION", "NAME", "STATE", "START_TIME", "END_TIME", "NODELIST"]
############### TODO : These should all be created from a configuration file at some point ###############


class PoolManager:
    """Class representing a single Run"""

    def __init__(self, **kwargs):
        self.site = config.conf.get("SITE", "id"),

    def add_worker_pool():
        logger.info("add_worker_pool")
        return None

    def rm_worker_pool():
        return None

    def get_current_slurm_workers(site: str = COMPUTE_SITE, ret_df: bool = False) -> Dict:
        slurm_status = {}

        _stdout, _stderr, _errcode = run_sh_command(condor_q_cmd, show_stdout=False)

        # NOTE: Leaving in for furute SuperFacility calls at cori
        # _stdout = sfapi.custom_cmd(
        #     token=access_token.token, cmd=squeue_cmd, site=site)
        # try:
        #     _stdout = _stdout['output']
        # except KeyError:
        #     logger.error("No output from squeue")
        #     output = {}
        #     for _type in compute_types:
        #         for _state in ["regular", "pending"]:
        #             output[f'{_type}_{_state}']
        #     return output

        # If we have an error return a dictionary with 0 for each type and state
        if _errcode != 0:
            logger.error(f"No output from squeue: {condor_q_cmd}")
            output = {}
            for _type in compute_types:
                for _state in ["regular", "pending"]:
                    output[f'{_type}_{_state}'] = 0
            return output

        # Gets jobs from output by splitting on new lines
        jobs = [job.split() for job in _stdout.split("\n")]

        # Make dataframe to query with
        df = pd.DataFrame(jobs, columns=squeue_columns)

        # replace N/A time with a value and convert times to datetime objects
        for time_col in ["START_TIME", "END_TIME"]:
            df.loc[df[time_col] == "N/A", time_col] = "2000-01-01T00:00:00"
            df[time_col] = pd.to_datetime(df[time_col])

        df["TIME_LEFT"] = df["END_TIME"] - df["START_TIME"]

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

        logger.info(slurm_status)
        return slurm_status

    # def get_condor_job_queue() -> pd.DataFrame:
    #     # Runs a condor_q autoformat to get the desired columns back
    #     _stdout, se, ec = run_sh_command(condor_q_cmd, show_stdout=False)
    #     if ec != 0:
    #         print(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
    #         return None

    #     # split outputs by rows
    #     outputs = _stdout.split("\n")

    #     # Get the column names from configuration
    #     columns = wanted_columns.split()
    #     # Split each row into columns
    #     queued_jobs = [job.split() for job in outputs]
    #     # Removes columns with no values (Usually the last column)
    #     queued_jobs = [q for q in queued_jobs if len(q) != 0]

    #     # Create a dataframe from the split outputs
    #     df = pd.DataFrame(queued_jobs, columns=columns)
    #     # Change the type
    #     df["RequestMemory"] = df["RequestMemory"].astype(int)
    #     df["JobStatus"] = df["JobStatus"].astype(int)
    #     df["RequestMemory"] = df["RequestMemory"].astype(float) / 1024
    #     df["RequestCpus"] = df["RequestCpus"].astype(float)
    #     df["CumulativeRemoteSysCpu"] = df["CumulativeRemoteSysCpu"].astype(float)
    #     df["CumulativeRemoteUserCpu"] = df["CumulativeRemoteUserCpu"].astype(float)

    #     now = int(time.time())
    #     df["JobStartDate"] = df["JobStartDate"].str.replace('undefined', str(now))
    #     df["total_running_time"] = now - df["JobStartDate"].astype(int)
    #     df["cpu_percentage"] = (((df['CumulativeRemoteSysCpu'] + df['CumulativeRemoteUserCpu']
    #                               ) / df['RequestCpus']) / df['total_running_time']) * 100

    #     df["total_q_time"] = df["JobStartDate"].astype(
    #         int) - df["QDate"].astype(int)

    #     return df

    # def determine_condor_job_sizes(df: pd.DataFrame):
    #     df["mem_bin"] = pd.cut(
    #         df["RequestMemory"],
    #         bins=[0, worker_sizes["regular_mem"], 10_000_000],
    #         labels=["regular", "over-mem"],
    #     )
    #     df["cpu_bin"] = pd.cut(
    #         df["RequestCpus"],
    #         bins=[0, worker_sizes["regular_cpu"], 10_000_000],
    #         labels=["regular", "over-cpu"],
    #     )

    #     condor_q_status = {}
    #     mask_running_status = df["JobStatus"].astype(int) == 2
    #     mask_idle_status = df["JobStatus"].astype(int) == 1
    #     mask_mem_regular = df["mem_bin"].str.contains("regular")
    #     mask_cpu_regular = df["cpu_bin"].str.contains("regular")

    #     mask_over = df["cpu_bin"].str.contains(
    #         "over") | df["mem_bin"].str.contains("over")

    #     mask_idle_regular = mask_mem_regular & mask_cpu_regular & mask_idle_status

    #     condor_q_status["idle_regular"] = sum(mask_idle_regular)
    #     condor_q_status["running_regular"] = sum(
    #         mask_mem_regular & mask_cpu_regular & mask_running_status
    #     )

    #     condor_q_status["hold_and_impossible"] = sum(mask_over)

    #     condor_q_status["regular_cpu_needed"] = sum(
    #         df[mask_idle_regular].RequestCpus)
    #     condor_q_status["regular_mem_needed"] = sum(
    #         df[mask_idle_regular].RequestMemory)

    #     return condor_q_status

    # def need_new_nodes(condor_job_queue: Dict, slurm_workers: Dict, machine: Dict) -> Dict:
    #     """
    #     Using the two dictionaries from the condor_q and squeue determine if we need any new workers fpr the machine type.
    #     """
    #     workers_needed = 0

    #     # Determines how many full (or partially full nodes) we need to create
    #     _cpu = (
    #         condor_job_queue[f"{machine}_cpu_needed"] /
    #         worker_sizes[f"{machine}_cpu"]
    #     )
    #     _mem = (
    #         condor_job_queue[f"{machine}_mem_needed"] /
    #         worker_sizes[f"{machine}_mem"]
    #     )
    #     _cpu = math.floor(_cpu)
    #     _mem = math.floor(_mem)

    #     workers_needed += max(_cpu, _mem)

    #     # If full workers_needed is 0 but we have work to be done still get a node
    #     if workers_needed == 0:
    #         if condor_job_queue[f"{machine}_cpu_needed"] or condor_job_queue[f"{machine}_mem_needed"]:
    #             workers_needed = 1

    #     # Total number running and pending to run (i.e. worker pool)
    #     current_pool_size = (
    #         slurm_workers[f"{machine}_pending"]
    #         + slurm_workers[f"{machine}_running"]
    #     )

    #     # If workers_needed is higher than the pool we'll add the diference
    #     # Else we don't need workers (add 0)
    #     workers_needed = max(0, workers_needed - current_pool_size)

    #     # If we have less running than the minimum we always need to add more
    #     # Either add what we need from queue (workers_needed)
    #     # Or what we're lacking in the pool (min - worker pool)
    #     if current_pool_size < MIN_POOL:
    #         workers_needed = max(MIN_POOL - current_pool_size, workers_needed)

    #     # Check to make sure we don't go over the max pool size
    #     if (workers_needed + current_pool_size) > MAX_POOL:
    #         # Only add up to max pool and no more
    #         workers_needed = MAX_POOL - current_pool_size

    #     return workers_needed

    # def auto_worker(site):
    #     if site == 'cori':
    #         ret = sfapi.post_job(token=access_token.token,
    #                              site='cori', script='/global/homes/t/tylern/spin_condor/worker.cori.job', isPath=True)
    #     elif site == 'perlmutter':
    #         ret = sfapi.post_job(token=access_token.token,
    #                              site=site, script='/global/homes/t/tylern/alvarez_condor/worker.perlmutter.ss11.job', isPath=True)
    #     else:
    #         return "Not a Valid Site"
    #     try:
    #         return repr(ret['jobid'])
    #     except Exception as e:
    #         return repr(e)

    # workers_needed = {}

    # def run_workers_needed():
    #     global job_queue_df
    #     job_queue_df = get_condor_job_queue()
    #     condor_job_queue = determine_condor_job_sizes(job_queue_df)
    #     logger.info(condor_job_queue)
    #     slurm_workers = get_current_slurm_workers(COMPUTE_SITE)
    #     logger.info(slurm_workers)
    #     global workers_needed
    #     workers_needed = {
    #         "regular": need_new_nodes(condor_job_queue, slurm_workers, "regular"),
    #     }

    #     logger.info(workers_needed)
    #     if workers_needed['regular'] > 0:
    #         jobid = auto_worker(COMPUTE_SITE)
    #         logger.info(f'Starting job {jobid}')

    #     return workers_needed

    # def run_cleanup():
    #     # Runs a condor_q autoformat to get the desired columns back
    #     _stdout, se, ec = run_sh_command(condor_idle_nodes, show_stdout=False)
    #     if ec != 0:
    #         print(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
    #         return None
    #     nodes = _stdout.split("\n")[:-1]
    #     logger.info(nodes)

    #     if len(nodes) == 0:
    #         return None
    #     try:
    #         logger.info(slurm_running_df.shape)
    #     except NameError as e:
    #         logger.info('no slurm nodes yet')
    #         return None

    #     for node in nodes:
    #         try:
    #             job_id = slurm_running_df[slurm_running_df['NODELIST(REASON)'] ==
    #                                       node].JOBID.iloc[0]
    #         except IndexError:
    #             continue

    #         logger.info(f"Removing {node} with JobID {job_id}")
    #         x = sfapi.delete_job(access_token.token,
    #                              site=COMPUTE_SITE, jobid=job_id)
    #         logger.info(x)

    #     return nodes


if __name__ == '__main__':
    pool = PoolManager()
    pool.get_current_slurm_workers()
