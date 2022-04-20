#!/usr/bin/env python
"""
Condor pool manager

Author: Seung-Jin Sul (ssul@lbl.gov)

Steps
1. Collect all RUNNING and PENDING SLURM job ids for Condor pool
2. Get each number of RUNNING and PENDING jobs per the memory requirement (normal, jgi_shared, jgi_exvivo)
3. Check if any of running SLURM jobs are using the node(s) and `scancel` them if not and the number of R and PD SLURM
   jobs is bigger than MIN pool size

"""
import click
import logging
import os
import time
import shlex
from datetime import datetime
import configparser
import jaws_condor.config
from jaws_condor import config
form jaws_condor.utils import run_sh_command


logger = None


def start_file_logger(filename, name='jaws_condor_pool_remove.log', level=logging.DEBUG, format_string=None):
    """Add a stream log handler.

    Args:
        - filename (string): Name of the file to write logs to
        - name (string): Logger name
        - level (logging.LEVEL): Set the logging level.
        - format_string (string): Set the format string

    Returns:
       -  None
    """
    if format_string is None:
        format_string = "%(asctime)s.%(msecs)03d %(name)s:%(lineno)d [%(levelname)s]  %(message)s"

    global logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string, datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def run_scancel(
    sq_cmd: str, running_condor_jobs: list, sc_cmd: str, keep_min_pool=False
):
    if len(running_condor_jobs) == 0:
        num_running_slurm_ids = []
        so, se, ec = run_sh_command(sq_cmd, show_stdout=False)
        if ec != 0:
            print(f"ERROR: failed to execute squeue command: {sq_cmd}")
            exit(1)
        num_running_slurm_ids = so.rstrip().split("\n")
        # Remove empty items
        num_running_slurm_ids = list(filter(None, num_running_slurm_ids))
        print(f"Collected SLURM jobs: {num_running_slurm_ids}")
        if keep_min_pool:
            # Keep min_pool_size number of nodes
            num_running_slurm_ids.sort()
            num_running_slurm_ids = num_running_slurm_ids[:-min_pool_size]
        print(sq_cmd)
        print(f"Candidate SLURM jobs to remove: {num_running_slurm_ids}")
        print("Number of SLURM jobs to remove: %d" % len(num_running_slurm_ids))
        if len(num_running_slurm_ids):
            sc_cmd = sc_cmd % ",".join(num_running_slurm_ids)
            so, se, ec = run_sh_command(sc_cmd, show_stdout=False)
            if ec != 0:
                print(f"ERROR: failed to execute scancel command: {sc_cmd}")
                exit(1)
            print(sc_cmd)
    else:
        print(f"Nothing to scancel from {sq_cmd}")


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=None,
                        help="Config file")
    parser.add_argument("-l", "--logfile", default=None,
                        help="Path to logfile")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")

    args = parser.parse_args()
    jaws_condor.config.Configuration(args.config)
    site_id = jaws_condor.config.conf.get_site_id()
    site_config = configparser.ConfigParser()
    if site_id == "CORI":
        site_config.read("condor_config/cori_config.ini")
    elif site_id == "JGI":
        site_config.read("condor_config/jgi_config.ini")
    elif site_id == "TAHOMA":
        site_config.read("condor_config/tahoma_config.ini")
    else:
        raise ValueError("Unknown side_id specified in condor backend config")

    #
    # Run condor_q_cmd to get the jobs in IDLE status
    #
    so, se, ec = run_sh_command(condor_q_cmd, show_stdout=False)
    if ec != 0:
        print(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
        exit(1)
    print("RUNNING Condor jobs")
    print("Job_id\tReq_mem\tReq_disk\tReq_cpu")
    print(f"{so.rstrip()}")

    normal_job = []
    shared_job = []
    exvivo_job = []

    for l in so.split("\n"):
        if l and len(shlex.split(l)) == 6:
            tok = shlex.split(l)
            job_id = tok[0]
            req_mem = float(tok[1])
            mem_unit = tok[2].strip()
            if mem_unit in ("KB", "kb"):
                req_mem = req_mem / (1024 * 1024)
            if mem_unit in ("MB", "mb"):
                req_mem = req_mem / 1024
            req_disk = float(tok[3])
            req_cpu = int(tok[5])

            if req_mem <= 118.0:
                normal_job.append([job_id, req_mem, req_disk, req_cpu])
            elif req_mem <= 740.0 and req_cpu <= 18:
                shared_job.append([job_id, req_mem, req_disk, req_cpu])
            elif req_mem <= 1450.0 and req_cpu <= 36:
                exvivo_job.append([job_id, req_mem, req_disk, req_cpu])
            else:
                print(
                    f"ERROR: no machines available for the job id, {job_id}, requested_memory={req_mem}, requested_disk={req_disk}, requested_cpu={req_cpu}"
                )

    print(f"Normal job: {normal_job}")
    print(f"Jgi_shared job: {shared_job}")
    print(f"Jgi_exvivo job: {exvivo_job}")

    #
    # Run scancel
    #
    run_scancel(squeue_cmd, normal_job, scancel_cmd, keep_min_pool=True)
    run_scancel(squeue_cmd_esslurm_shared, shared_job, scancel_cmd_esslurm)
    run_scancel(squeue_cmd_esslurm_exvivo, exvivo_job, scancel_cmd_esslurm)
