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
import argparse
import click
import logging
import os
import time
import shlex
from datetime import datetime
import configparser as cparser
import json
import jaws_condor.config
from jaws_condor import config
from jaws_condor.utils import run_sh_command


logger = None


def start_file_logger(
    filename,
    name="jaws_condor_pool_remove.log",
    level=logging.DEBUG,
    format_string=None,
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
    streamLogger = logging.StreamHandler()
    streamLogger.setLevel(level)
    streamLogger.setFormatter(formatter)
    logger.addHandler(streamLogger)


def run_scancel(
    rsc_t: str,
    sq_cmd: str,
    sc_cmd: str,
    running_condor_jobs: list,
    min_pool_sz: int,
    site_id=None,
):
    sq_cmd = sq_cmd.replace("<poolsz>", rsc_t) + " | awk '{ print $1; }'"
    if site_id is not None and site_id == "CORI" and rsc_t in ("large", "xlarge"):
        sq_cmd = "module load esslurm && " + sq_cmd
        sc_cmd = "module load esslurm && " + sc_cmd
    logger.debug(sq_cmd)
    logger.debug(sc_cmd)

    if len(running_condor_jobs) == 0:
        num_running_slurm_ids = []
        so, se, ec = run_sh_command(sq_cmd, log=logger, show_stdout=False)
        if ec != 0:
            logger.info(f"ERROR: failed to execute squeue command: {sq_cmd}")
            exit(1)
        num_running_slurm_ids = so.rstrip().split("\n")
        # Remove empty items
        num_running_slurm_ids = list(filter(None, num_running_slurm_ids))
        logger.info(f"Collected SLURM jobs: {num_running_slurm_ids}")
        if min_pool_sz > 0:
            # Keep min_pool_size number of nodes
            num_running_slurm_ids.sort()
            num_running_slurm_ids = num_running_slurm_ids[:-min_pool_sz]
        logger.info(sq_cmd)
        logger.info(f"Candidate SLURM jobs to remove: {num_running_slurm_ids}")
        logger.info("Number of SLURM jobs to remove: %d" % len(num_running_slurm_ids))
        if len(num_running_slurm_ids):
            sc_cmd = sc_cmd % ",".join(num_running_slurm_ids)
            so, se, ec = run_sh_command(sc_cmd, log=logger, show_stdout=False)
            if ec != 0:
                logger.info(f"ERROR: failed to execute scancel command: {sc_cmd}")
                exit(1)
            logger.info(sc_cmd)
    else:
        logger.info(f"Nothing to scancel from {sq_cmd}")


def collect_condor_running_jobs(condor_q_out: str, ram_range: list) -> dict:
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
    parser.add_argument("-s", "--site_config", default=None, help="Site config file")
    parser.add_argument("-l", "--logfile", default=None, help="Path to logfile")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enables debug logging"
    )

    args = parser.parse_args()
    jaws_condor.config.Configuration(args.config)
    site_id = jaws_condor.config.conf.get_site_id().upper()

    if args.logfile:
        logfile_path = logfile
    else:
        logfile_path = "./jaws_condor_remove.log"
    os.makedirs(os.path.dirname(logfile_path), exist_ok=True)
    start_file_logger(logfile_path, level=logging.DEBUG if args.debug else logging.INFO)

    site_config = cparser.ConfigParser(interpolation=cparser.ExtendedInterpolation())
    if site_id == "CORI":
        if args.site_config:
            site_config_path = args.site_config
        else:
            site_config_path = "cori_config.ini"
    elif site_id == "JGI":
        if args.site_config:
            site_config_path = args.site_config
        else:
            site_config_path = "jgi_config.ini"
    elif site_id == "TAHOMA":
        if args.site_config:
            site_config_path = args.site_config
        else:
            site_config_path = "tahoma_config.ini"
    else:
        raise ValueError("Unknown side_id specified in condor backend config")

    # Run condor_q_cmd to get the jobs in IDLE status
    site_config.read(site_config_path)
    condor_q_cmd = site_config["CONDOR"]["condor_q_cmd_rm"]
    so, se, ec = run_sh_command(condor_q_cmd, log=logger, show_stdout=False)
    if ec != 0:
        logger.critical(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
        exit(1)
    logger.info("RUNNING Condor jobs")
    logger.info("Job_id\tReq_mem\tReq_disk\tReq_cpu")
    logger.info(f"{so.rstrip()}")
    ram_range = json.loads(site_config.get("RESOURCE", "ram"))
    run_list = collect_condor_running_jobs(so, ram_range)
    logger.info(f"Running Condor job list: {run_list}")

    min_pool_size = json.loads(site_config.get("CONDOR", "min_pool"))
    resource_types = json.loads(site_config.get("RESOURCE", "types"))
    squeue_cmd_r_p_rm = site_config["SLURM"]["squeue_cmd_running_pending_rm"]
    scancel_cmd = site_config["SLURM"]["scancel_cmd"]

    # Run scancel
    for idx, t in enumerate(resource_types):
        logger.info(f"=== Checking for {t} ===")
        logger.debug(f"MIN pool size = {min_pool_size[idx]}")
        run_scancel(
            t,
            squeue_cmd_r_p_rm,
            scancel_cmd,
            run_list[idx],
            min_pool_size[idx],
            site_id=site_id,
        )
