#!/usr/bin/env python
"""
Condor pool manager for CORI
- condor_pool_remove.py

Author: Seung-Jin Sul (ssul@lbl.gov)

Steps
1. Collect all RUNNING and PENDING SLURM job ids for Condor pool
2. Get each number of RUNNING and PENDING jobs per the memeory requirement (normal, jgi_shared, jgi_exvivo)
3. Check if any of running SLURM jobs are using the node(s) and `scancel` them if not

"""
import os
import time 
import subprocess
import shlex
from utils import run_sh_command


#
# System vars, dirs, and cmds
#
accnt_name = "jaws_jtm"
condor_root = "/global/cfs/cdirs/jaws/condor"
condor_q_cmd = f"condor_q -pr {condor_root}/fmt_nobatch_running.cpf"
sbatch_cmd = "sbatch --parsable %s"
sbatch_cmd_esslurm = "module load esslurm && sbatch --parsable %s"
squeue_cmd = f"""squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t R,PD | grep jaws_condor_normal_worker """ + "| awk '{ print $1; }'"
squeue_cmd_esslurm_shared = f"""module load esslurm && squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t R,PD | grep jaws_condor_highmem_jgi_shared_worker""" + "| awk '{ print $1; }'"
squeue_cmd_esslurm_exvivo = f"""module load esslurm && squeue --format="%.18i %.9P %.40j %.8u %.10T %.10M %.9l %.6D %R" --me -u {accnt_name} -t R,PD | grep jaws_condor_highmem_jgi_exvivo_worker""" + "| awk '{ print $1; }'"
scancel_cmd = "scancel %s"
scancel_cmd_esslurm = "module load esslurm && scancel %s"


#
# Execute scancel
#
def run_scancel(sq_cmd: str, running_condor_jobs: list, scancel_cmd: str):
    if len(running_condor_jobs) == 0:
        num_running_slurm_ids = []
        so, se, ec = run_sh_command(sq_cmd, show_stdout=False)
        if ec != 0:
            print(f"ERROR: failed to execute squeue command: {sq_cmd}")
            exit(1)
        num_running_slurm_ids = so.rstrip().split('\n')
        num_running_slurm_ids = list(filter(None, num_running_slurm_ids))
        print(sq_cmd)
        print(f"RUNNING SLURM jobs: {num_running_slurm_ids}")
        print("Number of SLURM jobs to remove: %d" % len(num_running_slurm_ids))
        if len(num_running_slurm_ids):
            sc_cmd = scancel_cmd % ",".join(num_running_slurm_ids)
            so, se, ec = run_sh_command(sc_cmd, show_stdout=False)
            if ec != 0:
                print(f"ERROR: failed to execute scancel command: {sc_cmd}")
                exit(1)
            print(sc_cmd)
    else:
        print(f"Nothing to scancel from {sq_cmd}")


#
# Create condor print format file
#
with open(f"{condor_root}/fmt_nobatch_running.cpf", 'w') as cpf:
    cpf.write("""# condor_q format to get RUNNING and COMPLETED jobs
SELECT NOHEADER NOTITLE
   ClusterId     AS "    ID"  NOSUFFIX WIDTH AUTO
   RequestMemory AS REQUEST_MEMORY    WIDTH 10    PRINTAS READABLE_MB
   RequestDisk   AS REQUEST_DISK  WIDTH 12    PRINTAS READABLE_KB
   RequestCpus   AS REQUEST_CPUS
WHERE (JobStatus == 2 || JobStatus == 4)""")


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
            req_mem = req_mem/(1024*1024)
        if mem_unit in ("MB", "mb"):
            req_mem = req_mem/1024
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
run_scancel(squeue_cmd, normal_job, scancel_cmd)
run_scancel(squeue_cmd_esslurm_shared, shared_job, scancel_cmd_esslurm)
run_scancel(squeue_cmd_esslurm_exvivo, exvivo_job, scancel_cmd_esslurm)
