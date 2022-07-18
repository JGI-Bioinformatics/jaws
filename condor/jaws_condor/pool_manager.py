import logging
import time
import json
from jaws_site import config
import shlex
from jaws_condor.cmd_utils import run_sh_command, run_slurm_cmd, mem_unit_to_g

logger = logging.getLogger(__package__)


class PoolManager:
    """Class representing a single Run"""

    def __init__(self, **kwargs):
        self.config = config.conf.get_section("POOL_MANAGER")  # TODO Get right section
        self.site = self.config["SITE_ID"]

    def add_workers(self):
        logger.info("Checking to add workers to pool")

        condor_q_cmd = self.config["CONDOR"]["condor_q_cmd"]
        so, se, ec = run_sh_command(condor_q_cmd, log=logger, show_stdout=False)
        if ec != 0:
            logger.critical(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
            exit(1)
        logger.info("IDLE Condor jobs")
        logger.info("Job_id\tReq_mem\tReq_disk\tReq_cpu")
        logger.info(f"{so.rstrip()}")
        ram_range = json.loads(self.config.get("RESOURCE", "mem"))
        idle_list = collect_condor_jobs(so, ram_range)
        logger.info(f"IDLE Condor job list: {idle_list}")

        sbatch_cmd = self.config["SLURM"]["sbatch_cmd"]
        squeue_cmd_p = self.config["SLURM"]["squeue_cmd_pending"]
        squeue_cmd_r_p = self.config["SLURM"]["squeue_cmd_running_pending"]
        max_pool_size = json.loads(self.config.get("CONDOR", "max_pool"))
        min_pool_size = json.loads(self.config.get("CONDOR", "min_pool"))
        cpu_size = json.loads(self.config.get("RESOURCE", "cpu"))
        resource_types = json.loads(self.config.get("RESOURCE", "types"))
        job_files = [v for (k, v) in self.config.items("WORKER_TYPES")]

        # Run sbatch
        for idx, t in enumerate(resource_types):
            logger.info(f"=== Checking for {t} ===")
            run_sbatch(
                t,
                squeue_cmd_p,
                squeue_cmd_r_p,
                idle_list[idx],
                job_files[idx],
                sbatch_cmd,
                max_pool_size[idx],
                ram_range[idx],
                cpu_size[idx],
                site_id=self.site,
            )
            keep_min_pool(
                t,
                min_pool_size[idx],
                self.config["SLURM"]["squeue_cmd_running_pending"],
                job_files[idx],
                self.config["SLURM"]["sbatch_cmd"],
            )
        return None

    def rm_workers(self):
        logger.info("Checking to remove workers to pool")
        # Run condor_q_cmd to get the jobs in IDLE status
        self.config = config.conf.get_section("POOL_MANAGER")  # TODO Get right section

        condor_q_cmd = self.config["CONDOR"]["condor_q_cmd_rm"]
        so, se, ec = run_sh_command(condor_q_cmd, log=logger, show_stdout=False)
        if ec != 0:
            logger.critical(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
            exit(1)
        logger.info("RUNNING Condor jobs")
        logger.info("Job_id\tReq_mem\tReq_disk\tReq_cpu")
        logger.info(f"{so.rstrip()}")
        ram_range = json.loads(self.config.get("RESOURCE", "mem"))
        run_list = collect_condor_running_jobs(so, ram_range)
        logger.info(f"Running Condor job list: {run_list}")

        min_pool_size = json.loads(self.config.get("CONDOR", "min_pool"))
        resource_types = json.loads(self.config.get("RESOURCE", "types"))
        squeue_cmd_r_p_rm = self.config["SLURM"]["squeue_cmd_running_pending_rm"]
        scancel_cmd = self.config["SLURM"]["scancel_cmd"]

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
                site_id=self.site,
            )

        return None


def collect_condor_running_jobs(condor_q_out: str, ram_range: list) -> dict:
    idle_jobs = [[], [], [], []]  # small, med, large, xlarge
    for a_job in condor_q_out.split("\n"):
        if a_job and len(shlex.split(a_job)) == 6:
            tok = shlex.split(a_job)
            job_id = tok[0]
            req_mem = mem_unit_to_g(tok[2].strip(), float(tok[1]))
            req_disk = float(tok[3])
            req_cpu = int(tok[5])

            for idx, rr in enumerate(ram_range):
                ram_start = int(rr.split("-")[0])
                ram_end = int(rr.split("-")[1])
                if ram_end != 0 and ram_start < req_mem <= ram_end:
                    idle_jobs[idx].append([job_id, req_mem, req_disk, req_cpu])

    return idle_jobs


def run_scancel(rsc_t: str,
                sq_cmd: str,
                sc_cmd: str,
                running_condor_jobs: list,
                min_pool_sz: int,
                site_id=None):

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


def calculate_node_needed(idle_list: list, ram_range: str, ncpu: int) -> int:
    # 0-118 --> max mem for this type = 118
    ram_e = int(ram_range.split("-")[1])
    summed = [sum(x) for x in zip(*idle_list)]
    print(summed)
    sum_ram = summed[1]
    sum_cpu = summed[3]
    return max(int(sum_ram / ram_e), int(sum_cpu / ncpu)) + 1


def run_sbatch(
    rsc_t: str,
    sq_cmd: str,
    sq_r_pd_cmd: str,
    idle_jobs: list,
    batch_script: str,
    sb_cmd: str,
    max_pool_sz: int,
    ram_s_e: str,
    cpu_s: int,
    site_id=None,
):
    logger.info("Number of IDLE condor jobs: %d" % len(idle_jobs))
    sq_cmd = sq_cmd.replace("<poolsz>", rsc_t)
    logger.debug(sq_cmd)
    sq_r_pd_cmd = sq_r_pd_cmd.replace("<poolsz>", rsc_t)
    logger.debug(sq_r_pd_cmd)
    if len(idle_jobs):
        if site_id is not None and site_id == "CORI" and rsc_t in ("large", "xlarge"):
            sq_cmd = "module load esslurm && " + sq_cmd
            sb_cmd = "module load esslurm && " + sb_cmd
        num_pending_jobs = 0
        num_r_pd_jobs = 0
        num_pending_jobs = run_slurm_cmd(sq_cmd)
        assert num_pending_jobs != -1
        logger.info(sq_cmd)

        logger.info(f"Number of PENDING slurm jobs: {num_pending_jobs}")
        num_sbatches_new = calculate_node_needed(idle_jobs, ram_s_e, cpu_s)
        logger.info(f"num_sbatches_new = {num_sbatches_new}")
        # num_sbatches = len(idle_jobs) - int(num_pending_jobs)
        num_sbatches = num_sbatches_new - int(num_pending_jobs)
        if max_pool_sz > 0:
            num_r_pd_jobs = run_slurm_cmd(sq_r_pd_cmd)
            if (num_r_pd_jobs + num_sbatches) > max_pool_sz:
                num_sbatches = max_pool_sz - num_r_pd_jobs
                logger.info(f"Current pool size = {num_r_pd_jobs}")
                logger.info(f"MAX pool size = {max_pool_sz}")
                logger.info(f"Adjusted number of nodes to add = {num_sbatches}")
        if num_sbatches > 0:
            for _ in range(num_sbatches):
                logger.info(run_slurm_cmd(sb_cmd % batch_script))
                logger.debug(sb_cmd % batch_script)
                time.sleep(0.5)


def keep_min_pool(
    rsc_t: str, min_pool_sz: int, sq_cmd: str, batch_script: str, sb_cmd: str
):
    sq_cmd = sq_cmd.replace("<poolsz>", rsc_t)
    logger.debug(sq_cmd)
    num_total_r_pd_jobs = 0
    num_total_r_pd_jobs = run_slurm_cmd(sq_cmd)
    assert num_total_r_pd_jobs != -1
    logger.info(
        f"Number of total RUNNING and PENDING slurm jobs: {num_total_r_pd_jobs}"
    )
    if num_total_r_pd_jobs <= min_pool_sz:
        to_add = min_pool_sz - num_total_r_pd_jobs
        logger.info(f"Need to add {to_add} nodes to keep MIN size pool for {rsc_t}")
        sb_cmd = sb_cmd % batch_script
        for _ in range(to_add):
            so, se, ec = run_sh_command(sb_cmd)
            if ec != 0:
                logger.critical(f"ERROR: failed to execute sbatch command: {sb_cmd}")
                exit(1)
            logger.info(sb_cmd)
            time.sleep(0.5)


def collect_condor_jobs(condor_q_out: str, ram_range: list) -> list:
    idle_jobs = [[], [], [], []]  # small, med, large, xlarge
    for a_job in condor_q_out.split("\n"):
        if a_job and len(shlex.split(a_job)) == 6:
            # ex)
            # 190        1.0 TB    1024.0 MB  1
            # 191     1024.0 KB    1024.0 MB  16
            # 192      390.6 GB    1024.0 MB  4
            # 193       97.7 GB    1024.0 MB  10
            tok = shlex.split(a_job)
            job_id = tok[0]
            req_mem = mem_unit_to_g(tok[2].strip(), float(tok[1]))
            req_disk = float(tok[3])
            req_cpu = int(tok[5])
            for idx, rr in enumerate(ram_range):
                ram_start = int(rr.split("-")[0])
                ram_end = int(rr.split("-")[1])
                if ram_end != 0 and ram_start < req_mem <= ram_end:
                    idle_jobs[idx].append([int(job_id), req_mem, req_disk, req_cpu])

    return idle_jobs
