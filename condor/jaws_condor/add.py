#!/usr/bin/env python
"""
Condor pool manager

Author: Seung-Jin Sul (ssul@lbl.gov)

Steps
1. Collect the IDLE condor job IDs
2. Get the number of nodes required per memory requirement (small, medium, large, xlarge)
3. Execute sbatch if the current SLURM jobs + number of sbatches needed < MAX pool size
4. Check if MIN pool size is being maintained, if not add more nodes up to MIN

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
from jaws_condor.utils import run_sh_command, run_slurm_cmd, read_json


logger = None


def start_file_logger(
    filename, name="jaws_condor_pool_add.log", level=logging.DEBUG, format_string=None
):
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
        format_string = (
            "%(asctime)s.%(msecs)03d %(name)s:%(lineno)d [%(levelname)s]  %(message)s"
        )

    global logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def run_sbatch(
    sq_cmd: str,
    sq_r_pd_cmd: str,
    idle_jobs: list,
    batch_script: str,
    sb_cmd: str,
    site_id=None,
):
    logger.info("Number of IDLE condor jobs: %d" % len(idle_jobs))
    if len(idle_jobs):
        if site_id is not None and site_id == "CORI":
            sq_cmd = "module load esslurm && " + sq_cmd
            sb_cmd = "module load esslurm && " + sb_cmd
        num_pending_jobs = 0
        num_r_pd_jobs = 0
        num_pending_jobs = run_slurm_cmd(sq_cmd)
        assert num_pending_jobs != -1
        logger.info(sq_cmd)

        logger.info(f"Number of PENDING slurm jobs: {num_pending_jobs}")
        num_sbatches = len(idle_jobs) - int(num_pending_jobs)
        if sq_r_pd_cmd is not None:
            num_r_pd_jobs = run_slurm_cmd(sq_r_pd_cmd)
            if (num_r_pd_jobs + num_sbatches) > max_pool_size:
                num_sbatches = max_pool_size - num_r_pd_jobs
                logger.info(f"Current pool size = {num_r_pd_jobs}")
                logger.info(f"MAX pool size = {max_pool_size}")
                logger.info(f"Adjusted number of nodes to add = {num_sbatches}")
        for _ in range(num_sbatches):
            logger.info(run_slurm_cmd(sb_cmd % batch_script))
            time.sleep(0.5)


def keep_min_pool(min_pool_size: int, sq_cmd: str, batch_script: str, sb_cmd: str):
    num_total_r_pd_jobs = 0
    num_total_r_pd_jobs = run_slurm_cmd(sq_cmd)
    assert num_total_r_pd_jobs != -1
    logger.info(
        f"Number of total RUNNING and PENDING slurm jobs: {num_total_r_pd_jobs}"
    )
    if num_total_r_pd_jobs <= min_pool_size:
        to_add = min_pool_size - num_total_r_pd_jobs
        logger.info(f"Add {to_add} nodes to Keep MIN size pool")
        sb_cmd = sb_cmd % batch_script
        for _ in range(to_add):
            so, se, ec = run_sh_command(sb_cmd)
            if ec != 0:
                logger.info(f"ERROR: failed to execute sbatch command: {sb_cmd}")
                exit(1)
            logger.info(sb_cmd)
            time.sleep(0.5)


def collect_condor_jobs(condor_q_out: str, ram_range: list) -> dict:
    idle_jobs = [[], [], [], []]  # small, med, large, xlarge
    for l in condor_q_out.split("\n"):
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

            for idx, rr in enumerate(ram_range):
                ram_start = int(rr.split("-")[0])
                ram_end = int(rr.split("-")[1])
                if ram_end != 0 and ram_start < req_mem <= ram_end:
                    idle_jobs[idx].append([job_id, req_mem, req_disk, req_cpu])

    return idle_jobs


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=None, help="Config file")
    parser.add_argument("-l", "--logfile", default=None, help="Path to logfile")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enables debug logging"
    )

    args = parser.parse_args()
    jaws_condor.config.Configuration(args.config)
    site_id = jaws_condor.config.conf.get_site_id()
    site_config = configparser.ConfigParser(
        interpolation=configparser.ExtendedInterpolation()
    )
    if site_id == "CORI":
        site_config.read("condor_config/cori_config.ini")
    elif site_id == "JGI":
        site_config.read("condor_config/jgi_config.ini")
    elif site_id == "TAHOMA":
        site_config.read("condor_config/tahoma_config.ini")
    else:
        raise ValueError("Unknown side_id specified in condor backend config")

    condor_q_cmd = site_config["CONDOR"]["condor_q_cmd"]
    so, se, ec = run_sh_command(condor_q_cmd, log=logger, show_stdout=False)
    if ec != 0:
        logger.info(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
        exit(1)
    logger.info("IDLE Condor jobs")
    logger.info("Job_id\tReq_mem\tReq_disk\tReq_cpu")
    logger.info(f"{so.rstrip()}")

    ram_range = json.loads(site_config.get("RESOURCE", "ram"))
    idle_list = collect_condor_jobs(so, ram_range)
    logger.info(f"IDLE Condor job list: {idle_list}")

    sbatch_cmd = site_config["SLURM"]["sbatch_cmd"]
    squeue_cmd_p = site_config["SLURM"]["squeue_cmd_pending"]
    squeue_cmd_r_p = site_config["SLURM"]["squeue_cmd_running_pending"]

    run_sbatch(
        squeue_cmd_p % "small",
        None,
        idle_list[0],
        site_config["SLURM"]["small_slurm_job"],
        sbatch_cmd,
    )
    run_sbatch(
        squeue_cmd_p % "medium",
        squeue_cmd_r_p % "medium",
        idle_list[1],
        site_config["SLURM"]["medium_slurm_job"],
        sbatch_cmd,
    )
    run_sbatch(
        squeue_cmd_p % "large",
        None,
        idle_list[2],
        site_config["SLURM"]["large_slurm_job"],
        sbatch_cmd,
        site_id=site_id,
    )
    run_sbatch(
        squeue_cmd_p % "xlarge",
        None,
        idle_list[3],
        site_config["SLURM"]["xlarge_slurm_job"],
        sbatch_cmd,
        site_id=site_id,
    )

    # Maintain MIN pool size for medium compute nodes
    keep_min_pool(
        site_config.getint("CONDOR", "min_pool"),
        site_config["SLURM"]["squeue_cmd_running_pending"],
        site_config["SLURM"]["medium_slurm_job"],
        site_config["SLURM"]["sbatch_cmd"],
    )
