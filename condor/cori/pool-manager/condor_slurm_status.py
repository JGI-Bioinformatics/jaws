#!/usr/bin/env python
"""
Condor pool manager for CORI
- condor_pool_add.py

Author: Seung-Jin Sul (ssul@lbl.gov)

# Steps
# 1. Collect the IDLE condor job IDs
# 2. Get the number of nodes required per memory requirement (normal, jgi_shared, jgi_exvivo)
#    normal: <= 118G
#    jgi_shared: 118 ~ 740G
#    jgi_exvivo: 740 ~ 1450G
# 3. Execute sbatch

# Steps
1. Get idle and running htcondor jobs
    a. Match running htcondor slot to slurm worker
    b. See how many nodes can run idle htcondor jobs
        b.1 ) If no nodes can run task, add to "needed_node_{regular,highmem}" state
        b.2 ) If nodes are availible, add to "waiting_on_slurm_job" state
2. Collect current SLURM workers, pending/running
3. re-run step 1 to see if jobs are draining from htcondor queue
3. See what is needed and submit sbatch command 


"""

import shlex
from typing import Dict
from utils import run_sh_command

import pandas as pd


#
# System vars, dirs, and cmds
#
accnt_name = "jaws_jtm"
# condor_root = "/global/cfs/cdirs/jaws/condor"
condor_root = ".."
condor_q_cmd = f"condor_q -pr {condor_root}/fmt_nobatch_id.cpf"
normal_worker_q = f"{condor_root}/condor_worker_normal.job"
highmem_worker_jgishared_q = f"{condor_root}/condor_worker_highmem_jgi_shared.job"
highmem_worker_jgiexvivo_q = f"{condor_root}/condor_worker_highmem_jgi_exvivo.job"
sbatch_cmd = "sbatch --parsable %s"
sbatch_cmd_esslurm = "module load esslurm && sbatch --parsable %s"


MIN_POOL = 3
MAX_POOL = 100

#
# Execute sbatch
#


def get_current_slurm_workers():
    squeue_cmd = f"""squeue --clusters=all --format="%.18i %.24P %.100j %.8u %.10T %S %e" -p genepool,genepool_shared,exvivo,exvivo_shared"""

    _stdout, _stderr, error = run_sh_command(squeue_cmd, show_stdout=False)
    if error != 0:
        print(f"ERROR: failed to execute squeue command: {squeue_cmd}")
        ## return None

    # Gets jobs from output
    jobs = [job.split() for job in _stdout.split("\n") if len(job.split()) > 3]
    # Get column names from list
    columns = jobs[0]
    num_cols = len(columns)
    # Remove all instances of column names from list of jobs
    jobs = [job for job in jobs if job != columns and len(job) == num_cols]

    # Make dataframe to query with
    df = pd.DataFrame(jobs, columns=columns)

    # replace time with value and convert times to datetime
    for time_col in ['START_TIME', 'END_TIME']:
        df.loc[df[time_col] == 'N/A', time_col] = '2000-01-01T00:00:00'
        df[time_col]= pd.to_datetime(df[time_col])

    jaws_jtm_status = {}

    mask_jtm = df['USER'] == 'jaws_jtm'
    mask_genepool = df['PARTITION'].str.contains('genepool')
    mask_pending = df['STATE'] == 'PENDING'

    jaws_jtm_status['jaws_geenepool_pending'] = sum(mask_jtm & mask_pending & mask_genepool)
    jaws_jtm_status['jaws_geenepool_running'] = sum(mask_jtm & ~mask_pending & mask_genepool)

    jaws_jtm_status['jaws_exvivo_pending'] = sum(mask_jtm & mask_pending & ~mask_genepool)
    jaws_jtm_status['jaws_exvivo_running'] = sum(mask_jtm & ~mask_pending & ~mask_genepool)

    jaws_jtm_status['other_genepool_pending'] = sum(~mask_jtm & mask_pending & mask_genepool)
    jaws_jtm_status['other_genepool_running'] = sum(~mask_jtm & ~mask_pending & mask_genepool)

    jaws_jtm_status['other_exvivo_pending'] = sum(~mask_jtm & mask_pending & ~mask_genepool)
    jaws_jtm_status['other_exvivo_running'] = sum(~mask_jtm & ~mask_pending & ~mask_genepool)

    return jaws_jtm_status

def get_needed_condor_workers() -> Dict:
    #
    # Create condor print format file
    #
    with open(f"{condor_root}/fmt_nobatch_id.cpf", 'w') as cpf:
        cpf.write("""# condor_q format to list IDLE jobs
SELECT BARE
    ClusterId     AS "    ID"  NOSUFFIX WIDTH AUTO
    RequestMemory AS REQUEST_MEMORY    WIDTH AUTO    PRINTAS READABLE_MB
    RequestCpus   AS REQUEST_CPUS
    JobStatus     AS ST
""")
    
    ## Job status 
    ## 1 is idle
    ## 2 is running
    ## 5 is hold


    #
    # Run condor_q_cmd to get the jobs in IDLE status
    #
    _stdout, se, ec = run_sh_command(condor_q_cmd, show_stdout=False)
    if ec != 0:
        print(f"ERROR: failed to execute condor_q command: {condor_q_cmd}")
        exit(1)

    outputs = _stdout.split("\n")
    columns = ['ClusterId','RequestMemory', 'mem_unit','RequestCpus','JobStatus']
    queued_jobs = [job.split() for job in outputs]
    queued_jobs = [q for q in queued_jobs if len(q) != 0]
    df = pd.DataFrame(queued_jobs, columns=columns)
    df['mem_bin'] = pd.cut(df['RequestMemory'].astype(float), bins=[0, 120, 1450, 10_000_000], labels=['regular', 'exvivo', 'over-mem'])
    df['cpu_bin'] = pd.cut(df['RequestCpus'].astype(float), bins=[0, 64, 72, 256, 10_000_000], labels=['regular', 'exvivo', 'perlmutter', 'over-cpu'])

    condor_q_status = {}
    mask_running_status = (df['JobStatus'].astype(int) == 2)
    mask_idle_status = (df['JobStatus'].astype(int) == 1)
    mask_mem_regular = df['mem_bin'].str.contains("regular")
    mask_mem_exvivo = df['mem_bin'].str.contains("exvivo")
    mask_cpu_regular = df['cpu_bin'].str.contains("regular")
    mask_cpu_exvivo = df['cpu_bin'].str.contains("exvivo")

    mask_over = (df['cpu_bin'].str.contains("over") | df['mem_bin'].str.contains("over"))
    condor_q_status['idle_regular'] = sum(mask_mem_regular & mask_cpu_regular & mask_idle_status)
    condor_q_status['running_regular'] = sum(mask_mem_regular & mask_cpu_regular & mask_running_status)
    condor_q_status['idle_exvivo'] = sum((mask_mem_exvivo | mask_cpu_exvivo) & mask_idle_status)
    condor_q_status['running_exvivo'] = sum((mask_mem_exvivo | mask_cpu_exvivo) & mask_running_status)
    condor_q_status['impossible'] = sum(mask_over)

    return condor_q_status


if __name__ == '__main__':
    slurm_workers = get_current_slurm_workers()
    print(slurm_workers)
    condor_workers = get_needed_condor_workers()
    print(condor_workers)
