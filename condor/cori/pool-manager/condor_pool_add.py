#!/usr/bin/env python
"""
Condor pool manager for CORI
- condor_pool_add.py

Author: Seung-Jin Sul (ssul@lbl.gov)

Steps
1. Collect the IDLE condor job IDs
2. Get the number of nodes required per memory requirement (normal, jgi_shared, jgi_exvivo)
   normal: <= 118G
   jgi_shared: 118 ~ 740G
   jgi_exvivo: 740 ~ 1450G
3. Execute sbatch if the current SLURM jobs + number of sbatches needed < MAX pool size
4. Check if MIN pool size is being maintained, if not add more nodes up to MIN

"""
import os
import time
import subprocess
import shlex
from utils import run_sh_command, run_slurm_cmd
import configparser
from datetime import datetime

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))


#
# System vars, dirs, and cmds
#
config = configparser.ConfigParser()
config.read("jaws_condor_pool_manager.ini")
accnt_name = config["SLURM"]["accnt_name"]
condor_root = config["CONDOR"]["condor_root"]
min_pool_size = config.getint("CONDOR", "min_pool_size")
max_pool_size = config.getint("CONDOR", "max_pool_size")
condor_q_cmd = f"condor_q -pr {condor_root}/fmt_nobatch_id.cpf"
normal_worker_job = f"{condor_root}/condor_worker_normal.job"
highmem_worker_jgishared_job = f"{condor_root}/condor_worker_highmem_jgi_shared.job"
highmem_worker_jgiexvivo_job = f"{condor_root}/condor_worker_highmem_jgi_exvivo.job"
sbatch_cmd = "sbatch --parsable %s"
sbatch_cmd_esslurm = "module load esslurm && sbatch --parsable %s"
squeue_cmd = f"""squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t PENDING | grep jaws_condor_normal_worker | wc -l"""
squeue_cmd_running_pending = f"""squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t R,PD | grep jaws_condor_normal_worker | wc -l"""
squeue_cmd_esslurm_shared = f"""module load esslurm && squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t PENDING | grep jaws_condor_highmem_jgi_shared_worker | wc -l"""
squeue_cmd_esslurm_exvivo = f"""module load esslurm && squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t PENDING | grep jaws_condor_highmem_jgi_exvivo_worker | wc -l"""


#
# Create condor print format file
#
with open(f"{condor_root}/fmt_nobatch_id.cpf", "w") as cpf:
    cpf.write(
        """# condor_q format to list IDLE jobs
SELECT NOHEADER NOTITLE
   ClusterId     AS "    ID"  NOSUFFIX WIDTH AUTO
   RequestMemory AS REQUEST_MEMORY    WIDTH 10    PRINTAS READABLE_MB
   RequestDisk   AS REQUEST_DISK  WIDTH 12    PRINTAS READABLE_KB
   RequestCpus   AS REQUEST_CPUS
WHERE JobStatus == 1"""
    )


#
# Execute sbatch
#
def run_sbatch(
    sq_cmd: str, sq_r_pd_cmd: str, idle_jobs: list, batch_script: str, sb_cmd: str
):
    num_pending_jobs = 0
    num_r_pd_jobs = 0
    num_pending_jobs = run_slurm_cmd(sq_cmd)
    assert num_pending_jobs != -1
    print(sq_cmd)
    print("Number of IDLE condor jobs: %d" % len(idle_jobs))
    print(f"Number of PENDING slurm jobs: {num_pending_jobs}")
    num_sbatches = len(idle_jobs) - int(num_pending_jobs)
    if sq_r_pd_cmd is not None:
        num_r_pd_jobs = run_slurm_cmd(sq_r_pd_cmd)
        if (num_r_pd_jobs + num_sbatches) > max_pool_size:
            num_sbatches = max_pool_size - num_r_pd_jobs
            print(f"Current pool size = {num_r_pd_jobs}")
            print(f"MAX pool size = {max_pool_size}")
            print(f"Adjusted number of nodes to add = {num_sbatches}")
    for _ in range(num_sbatches):
        print(run_slurm_cmd(sb_cmd % batch_script))
        time.sleep(0.5)


#
# Check MIN pool
#
def keep_min_pool(sq_cmd: str, batch_script: str, sb_cmd: str):
    num_total_r_pd_jobs = 0
    num_total_r_pd_jobs = run_slurm_cmd(sq_cmd)
    assert num_total_r_pd_jobs != -1
<<<<<<< HEAD
    print(f"Number of total RUNNING and PENDING slurm jobs: {num_total_r_pd_jobs}")
=======
>>>>>>> caeeaa1b001c339040bce943a4db192b530383c4
    if num_total_r_pd_jobs <= min_pool_size:
        to_add = min_pool_size - num_total_r_pd_jobs
        print(f"Add {to_add} nodes to Keep MIN size pool")
        sb_cmd = sb_cmd % batch_script
        for _ in range(to_add):
            so, se, ec = run_sh_command(sb_cmd)
            if ec != 0:
                print(f"ERROR: failed to execute sbatch command: {sb_cmd}")
                exit(1)
            print(sb_cmd)
            time.sleep(0.5)


#
# Run condor_q_cmd to get the jobs in IDLE status
#
so, se, ec = run_sh_command(condor_q_cmd, show_stdout=False)
if ec != 0:
    print(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
    exit(1)
print("IDLE Condor jobs")
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
# Run sbatch
#
run_sbatch(
    squeue_cmd, squeue_cmd_running_pending, normal_job, normal_worker_job, sbatch_cmd
)
run_sbatch(
    squeue_cmd_esslurm_shared,
    None,
    shared_job,
    highmem_worker_jgishared_job,
    sbatch_cmd_esslurm,
)
run_sbatch(
    squeue_cmd_esslurm_exvivo,
    None,
    exvivo_job,
    highmem_worker_jgiexvivo_job,
    sbatch_cmd_esslurm,
)
keep_min_pool(squeue_cmd_running_pending, normal_worker_job, sbatch_cmd)
