#!/usr/bin/env python
"""

# Steps
1. Get idle and running htcondor jobs
    a. Match running htcondor slot to slurm worker
    b. See how many nodes can run idle htcondor jobs
        b.1 ) If no nodes can run task, add to "needed_node_{regular,highmem}" state
        b.2 ) If nodes are availible, add to "waiting_on_slurm_job" state
2. Collect current SLURM workers, pending/running
3. See what is needed and submit sbatch command 

"""
from typing import Dict
from utils import run_sh_command
import pandas as pd
import time
from pathlib import Path
import argparse


try:
    from SuperfacilityAPI import SuperfacilityAPI, SuperfacilityAccessToken
    superfacility = Path.home() / '.superfacility'
    pem = list(superfacility.glob('*.pem'))
    if len(pem) > 0:
        access_token = SuperfacilityAccessToken(key_path=pem[0])
        sfapi = SuperfacilityAPI(token=access_token.token)
except Exception as e:
    print(f"{type(e).__name__} : {e}")
    sfapi = None

###################### TODO : These should all be created from a configuration file at some point ######################
accnt_name = "jaws_jtm"
wanted_columns = "ClusterId RequestMemory RequestCpus JobStatus NumRestarts QDate"
condor_q_cmd = f"condor_q -af {wanted_columns}"

MIN_POOL = 3
MAX_POOL = 100

worker_sizes = {
    "regular_cpu": 64,
    "regular_mem": 128,
    "large_cpu": 72,
    "large_mem": 1450,
    "perlmutter_cpu": 256,
    "perlmutter_mem": 512,
}

default_form = {}
default_form["regular"] = {
    "site": "cori",
    "qos": "genepool_special",
    "constraint": "haswell",
    "account": "fungalp",
    "time": "06:00:00",
    "cluster": "cori",
    "script_location": ".",
}
default_form["large"] = {
    "site": "cori",
    "qos": "exvivo",
    "constraint": "skylake",
    "account": "fungalp",
    "time": "06:00:00",
    "cluster": "escori",
    "script_location": ".",
}

# Extra partition arguments for cori
squeue_args = {}
squeue_args["cori"] = "--clusters=all -p genepool,genepool_shared,exvivo,exvivo_shared"
squeue_args["jgi"] = "-p lr3"
squeue_args["tahoma"] = ""
############### TODO : These should all be created from a configuration file at some point ###############


def get_condor_job_queue() -> pd.DataFrame:
    # Runs a condor_q autoformat to get the desired columns back
    _stdout, se, ec = run_sh_command(condor_q_cmd, show_stdout=False)
    if ec != 0:
        print(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
        return None

    # split outputs by rows
    outputs = _stdout.split("\n")

    # Get the column names from configuration
    columns = wanted_columns.split()
    # Split each row into columns
    queued_jobs = [job.split() for job in outputs]
    # Removes columns with no values (Usually the last column)
    queued_jobs = [q for q in queued_jobs if len(q) != 0]

    # Create a dataframe from the split outputs
    df = pd.DataFrame(queued_jobs, columns=columns)
    # Change the type 
    df["RequestMemory"] = df["RequestMemory"].astype(int)
    df["JobStatus"] = df["JobStatus"].astype(int)
    df["RequestMemory"] = df["RequestMemory"].astype(float)
    df["RequestCpus"] = df["RequestCpus"].astype(float)
    df["total_q_time"] = int(time.time()) - df["QDate"].astype(int)

    return df

def determine_condor_job_sizes(df: pd.DataFrame):
    """
    Bins the jobs in htcondor queue into memory bins and cpu bins to determine which type of node should be used for each job.
    """

    # Use pd cut to bin from 0GB -> regular_memGB -> large_memGB -> Over memory
    df["mem_bin"] = pd.cut(
        df["RequestMemory"],
        bins=[0, worker_sizes["regular_mem"], worker_sizes["large_mem"], 10_000_000],
        labels=["regular", "large", "over-mem"],
    )
    # Use pd cut to bin from 0 cores -> regular_mem cores -> large_mem cores -> Over number of cores
    df["cpu_bin"] = pd.cut(
        df["RequestCpus"],
        bins=[0, worker_sizes["regular_cpu"], worker_sizes["large_cpu"], 10_000_000],
        labels=["regular", "large", "over-cpu"],
    )

    condor_q_status = {}
    # Makes masks based on the condor_q dataframe
    mask_running_status = df["JobStatus"].astype(int) == 2
    mask_idle_status = df["JobStatus"].astype(int) == 1
    mask_mem_regular = df["mem_bin"].str.contains("regular")
    mask_mem_large = df["mem_bin"].str.contains("large")
    mask_cpu_regular = df["cpu_bin"].str.contains("regular")
    mask_cpu_large = df["cpu_bin"].str.contains("large")

    # If either cpu/mem are over the limits we should do something with the job
    mask_over = df["cpu_bin"].str.contains("over") | df["mem_bin"].str.contains("over")

    mask_regular =  mask_mem_regular & mask_cpu_regular
    mask_large =  (mask_mem_large | mask_cpu_large)

    # If it sits in the regular bin for cpu and memorry
    mask_idle_regular = mask_regular & mask_idle_status
    # If either the memory or cpu need a larger node
    mask_idle_large = mask_large & mask_idle_status

    condor_q_status["idle_regular"] = sum(mask_regular & mask_idle_status)
    condor_q_status["running_regular"] = sum(mask_regular & mask_running_status)
    condor_q_status["idle_large"] = sum(mask_idle_large)
    condor_q_status["running_large"] = sum(
        mask_large & mask_running_status
    )
    condor_q_status["impossible"] = sum(mask_over)

    condor_q_status["regular_cpu_needed"] = sum(df[mask_idle_regular].RequestCpus)
    condor_q_status["regular_mem_needed"] = sum(df[mask_idle_regular].RequestMemory)

    condor_q_status["large_cpu_needed"] = sum(df[mask_idle_large].RequestCpus)
    condor_q_status["large_mem_needed"] = sum(df[mask_idle_large].RequestMemory)

    return condor_q_status


def get_current_slurm_workers_cmd(site: str = "cori") -> Dict:
    """
    Get outputs from slurm using a command line call

    Returns:
       Str
    """
    squeue_cmd = (
        f'squeue --format="%.18i %D %.24P %.100j %.20u %.10T %S %e" {squeue_args[site]}'
    )
    _stdout, _stderr, error = run_sh_command(squeue_cmd, show_stdout=False)
    if error != 0:
        print(f"ERROR: failed to execute squeue command: {squeue_cmd}")
        return None

    return _stdout

# returns the same as above but using sfapi
def get_current_slurm_workers_sfapi(site: str = "cori") -> Dict:
    """
    Get outputs from slurm using a sfapi call

    Returns:
       Str
    """
    squeue_cmd = f'squeue --format="%.18i %D %.24P %.100j %.20u %.10T %S %e" {squeue_args[site]}'
    try:
        _stdout = sfapi.custom_cmd(cmd=squeue_cmd)
        _stdout = _stdout['output']
        return _stdout
    except Exception as e:
        print(f"{type(e).__name__} : {e}")
        return None

def get_slurm_dataframe(site: str = "cori"):
    _stdout = get_current_slurm_workers_cmd()
    # _stdout = get_current_slurm_workers_sfapi()

    # Splits jobs from output
    jobs = [job.split() for job in _stdout.split("\n") if len(job.split()) > 3]

    # Get column names from list
    columns = jobs[0]
    num_cols = len(columns)
    # Remove all instances of column names from list of jobs
    jobs = [job for job in jobs if job != columns and len(job) == num_cols]

    # Convert the list into a dataframe
    df = pd.DataFrame(jobs, columns=columns)
    # Change types of the columns in the dataframe for use later
    df["NODES"] = df["NODES"].astype(int)

    # replace N/A time with value and convert times to datetime
    # Getting the start and end times could be useful later when
    # determining which nodes to remove from pool
    for time_col in ["START_TIME", "END_TIME"]:
        df.loc[df[time_col] == "N/A", time_col] = "2000-01-01T00:00:00"
        df[time_col] = pd.to_datetime(df[time_col])

    # Get's time left in the job
    df["TIME_LEFT"] = df["END_TIME"] - df["START_TIME"]

    return df

def get_current_slurm_workers(df):
    # Make a blank dictionary to return at the end
    jaws_jtm_status = {}

    # Create masks (list of true flase) depending on a selection criteria
    # We can combine these with logical and (&) or (|) not (~) to get desired selection
    mask_jaws = df["USER"] == accnt_name
    mask_regular = df["NAME"].str.contains("regular")
    mask_large = df["NAME"].str.contains("large")
    mask_pending = df["STATE"] == "PENDING"
    mask_running = df["STATE"] == "RUNNING"
    mask_condor = df["NAME"].str.contains("condor")

    # The jobs we care about are from jaws user with condor in the name
    mask_jaws = mask_jaws & mask_condor

    # Each of these selects for a certian type of node based on a set of masks
    # Add the number of nodes to get how many nodes are in each catogory
    jaws_jtm_status["jaws_regular_pending"] = sum(
        df[mask_jaws & mask_pending & mask_regular].NODES
    )
    jaws_jtm_status["jaws_regular_running"] = sum(
        df[mask_jaws & mask_running & mask_regular].NODES
    )

    jaws_jtm_status["jaws_large_pending"] = sum(
        df[mask_jaws & mask_pending & mask_large].NODES
    )
    jaws_jtm_status["jaws_large_running"] = sum(
        df[mask_jaws & ~mask_pending & mask_large].NODES
    )

    return jaws_jtm_status


def need_new_nodes(condor_job_queue: Dict, slurm_workers: Dict, machine: Dict) -> Dict:
    """
    Using the two dictionaries from the condor_q and squeue determine if we need any new workers fpr the machine type.
    """
    workers_needed = 0

    # Determines how many full (or partially full nodes) we need to create
    _cpu = (
        condor_job_queue[f"{machine}_cpu_needed"] // worker_sizes[f"{machine}_cpu"]
    )
    _mem = (
        condor_job_queue[f"{machine}_mem_needed"] // worker_sizes[f"{machine}_mem"]
    )

    workers_needed += max(_cpu, _mem)

    # If full workers_needed is 0 but we have work to be done still get a node
    if workers_needed == 0:
        if condor_job_queue[f"{machine}_cpu_needed"] or condor_job_queue[f"{machine}_mem_needed"]:
            workers_needed = 1

    # Total number running and pending to run (i.e. worker pool)
    current_pool_size = (
        slurm_workers[f"jaws_{machine}_pending"]
        + slurm_workers[f"jaws_{machine}_running"]
    )

    # If workers_needed is higher than the pool we'll add the diference
    # Else we don't need workers (add 0)
    workers_needed = max(0, workers_needed - current_pool_size)

    # If we have less running than the minimum we always need to add more
    # Either add what we need from queue (workers_needed)
    # Or what we're lacking in the pool (min - worker pool)
    if current_pool_size < MIN_POOL:
        workers_needed = max(MIN_POOL - current_pool_size, workers_needed)

    # Check to make sure we don't go over the max pool size
    if (workers_needed + current_pool_size) > MAX_POOL:
        # Only add up to max pool and no more
        workers_needed = MAX_POOL - current_pool_size

    return workers_needed

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=None, help="Config file")
    parser.add_argument("-l", "--logfile", default=None, help="Path to logfile")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enables debug logging"
    )

    args = parser.parse_args()

    print("============= condor_q =============")
    condor_dataframe = get_condor_job_queue()
     # condor_dataframe.to_html() -> Could be populated to a web frontend
    print(condor_dataframe)
    print("============= condor_q =============\n\n")
    print("============= condor_job_sizes =============")
    condor_job_sizes = determine_condor_job_sizes(condor_dataframe)
    print(condor_job_sizes)
    print("============= condor_job_sizes =============\n\n")
    print("============= slurm_jobs =============")
    slurm_dataframe = get_slurm_dataframe()
    # slurm_dataframe.to_html() -> Could be populated to a web frontend
    slurm_workers = get_current_slurm_workers(slurm_dataframe)
    print(slurm_workers)
    print("============= slurm_jobs =============\n\n")
    print("============= workers_needed =============")
    workers_needed = {_type : need_new_nodes(condor_job_sizes, slurm_workers, _type) for _type in ['regular', 'large']}
    print(workers_needed)
    print("============= workers_needed =============\n\n")


if __name__ == '__main__':
    cli()