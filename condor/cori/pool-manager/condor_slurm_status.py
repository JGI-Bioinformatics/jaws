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

import shlex
from typing import Dict
import uuid
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


###################### TODO : These should all be created from a configuration file at some point ######################
MIN_POOL = 3
MAX_POOL = 100

worker_sizes = {'regular_cpu' : 64, 'regular_mem':128, 'exvivo_cpu': 72, 'exvivo_mem':1450, 'perlmutter_cpu': 256, 'perlmutter_mem': 512}

default_form = {}
default_form['regular'] = {
                'qos': 'genepool_special',
                'constraint': 'haswell',
                'account': 'fungalp',
                'time': '06:00:00',
                'cluster': 'cori',
                'script_location': '.',
            }
default_form['exvivo'] = {
                'qos': 'exvivo',
                'constraint': 'skylake',
                'account': 'fungalp',
                'time': '06:00:00',
                'cluster': 'escori',
                'script_location': '.',
            }
############### TODO : These should all be created from a configuration file at some point ###############


def get_current_slurm_workers() -> Dict:
    partitions = '-p genepool,genepool_shared,exvivo,exvivo_shared'
    squeue_cmd = f'squeue --clusters=all --format="%.18i %D %.24P %.100j %.20u %.10T %S %e" {partitions}'

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
    df['NODES'] = df['NODES'].astype(int)

    # replace time with value and convert times to datetime
    for time_col in ['START_TIME', 'END_TIME']:
        df.loc[df[time_col] == 'N/A', time_col] = '2000-01-01T00:00:00'
        df[time_col]= pd.to_datetime(df[time_col])

    df["TIME_LEFT"] = df['END_TIME'] - df['START_TIME']
    jaws_jtm_status = {}

    mask_jtm = df['USER'] == 'jaws_jtm' 
    mask_genepool = df['PARTITION'].str.contains('genepool')
    mask_pending = df['STATE'] == 'PENDING'
    mask_condor = df['NAME'].str.contains('condor')

    mask_jtm = (mask_jtm & mask_condor)

    # Each of these selects for a certian type of node based on a set of masks
    # Add the number of nodes to get how many are in each catogory
    jaws_jtm_status['jaws_regular_pending'] = sum(df[mask_jtm & mask_pending & mask_genepool].NODES)
    jaws_jtm_status['jaws_regular_running'] = sum(df[mask_jtm & ~mask_pending & mask_genepool].NODES)

    jaws_jtm_status['jaws_exvivo_pending'] = sum(df[mask_jtm & mask_pending & ~mask_genepool].NODES)
    jaws_jtm_status['jaws_exvivo_running'] = sum(df[mask_jtm & ~mask_pending & ~mask_genepool].NODES)

    jaws_jtm_status['other_genepool_pending'] = sum(df[~mask_jtm & mask_pending & mask_genepool].NODES)
    jaws_jtm_status['other_genepool_running'] = sum(df[~mask_jtm & ~mask_pending & mask_genepool].NODES)

    jaws_jtm_status['other_exvivo_pending'] = sum(df[~mask_jtm & mask_pending & ~mask_genepool].NODES)
    jaws_jtm_status['other_exvivo_running'] = sum(df[~mask_jtm & ~mask_pending & ~mask_genepool].NODES)

    return jaws_jtm_status

def get_condor_job_queue() -> Dict:
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
    df['RequestMemory'] = df['RequestMemory'].astype(float)
    df['RequestCpus'] = df['RequestCpus'].astype(float)

    df['mem_bin'] = pd.cut(df['RequestMemory'], 
        bins=[0, worker_sizes['regular_mem'], worker_sizes['exvivo_mem'], 10_000_000], 
        labels=['regular', 'exvivo', 'over-mem'])
    df['cpu_bin'] = pd.cut(df['RequestCpus'], 
        bins=[0, worker_sizes['regular_cpu'], worker_sizes['exvivo_cpu'], 10_000_000], 
        labels=['regular', 'exvivo', 'over-cpu'])

    condor_q_status = {}
    mask_running_status = (df['JobStatus'].astype(int) == 2)
    mask_idle_status = (df['JobStatus'].astype(int) == 1)
    mask_mem_regular = df['mem_bin'].str.contains("regular")
    mask_mem_exvivo = df['mem_bin'].str.contains("exvivo")
    mask_cpu_regular = df['cpu_bin'].str.contains("regular")
    mask_cpu_exvivo = df['cpu_bin'].str.contains("exvivo")

    mask_over = (df['cpu_bin'].str.contains("over") | df['mem_bin'].str.contains("over"))

    mask_idle_regular = (mask_mem_regular & mask_cpu_regular & mask_idle_status)
    mask_idle_exvivo = ((mask_mem_exvivo | mask_cpu_exvivo) & mask_idle_status)

    condor_q_status['idle_regular'] = sum(mask_idle_regular)
    condor_q_status['running_regular'] = sum(mask_mem_regular & mask_cpu_regular & mask_running_status)
    condor_q_status['idle_exvivo'] = sum(mask_idle_exvivo)
    condor_q_status['running_exvivo'] = sum((mask_mem_exvivo | mask_cpu_exvivo) & mask_running_status)
    condor_q_status['hold_and_impossible'] = sum(mask_over)

    condor_q_status['regular_cpu_needed'] = sum(df[mask_idle_regular].RequestCpus)
    condor_q_status['regular_mem_needed'] = sum(df[mask_idle_regular].RequestMemory)

    condor_q_status['exvivo_cpu_needed'] = sum(df[mask_idle_exvivo].RequestCpus)
    condor_q_status['exvivo_mem_needed'] = sum(df[mask_idle_exvivo].RequestMemory)


    return condor_q_status


def need_new_nodes(condor_job_queue:Dict, slurm_workers:Dict, machine:Dict) -> Dict:
    workers_needed = 0
    
    # If we need more than a node add a node (or more)
    if condor_job_queue[f'{machine}_cpu_needed'] > worker_sizes[f'{machine}_cpu'] \
        or condor_job_queue[f'{machine}_mem_needed'] > worker_sizes[f'{machine}_mem']:
        # Calculate what we need
        _cpu = condor_job_queue[f'{machine}_cpu_needed']//worker_sizes[f'{machine}_cpu']
        _mem = condor_job_queue[f'{machine}_mem_needed']//worker_sizes[f'{machine}_mem']
        workers_needed += max(_cpu,_mem)

    # Total number running and pending to run (i.e. worker pool)
    current_pool_size = slurm_workers[f'jaws_{machine}_pending'] \
        + slurm_workers[f'jaws_{machine}_running']

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

job_template="""#!/bin/bash
#SBATCH -N {num_nodes}
#SBATCH -q {qos}
#SBATCH -C {constraint}
#SBATCH -A {account}
#SBATCH -t {time}
#SBATCH -M {cluster}
#SBATCH --job-name=jaws_condor_{machine}_worker
#SBATCH --exclusive
#SBATCH -e htcondor_%j.err
#SBATCH -o htcondor_%j.out

export WORKER_SCRIPTS="{script_location}"

# For each node start a worker on that node
for node in $(scontrol show hostnames ${{SLURM_NODELIST}}); do
    srun -N 1 -n 1 -c 1 --overlap ${{WORKER_SCRIPTS}}/start_worker.sh &
    sleep 2
done

# 6 hours is 21600 seconds
# Sleeps for a minute less then that to cleanup workers
# This is better for htcondor to handle instead of leaving hanging slots
# sleep 21540
sleep {time_sec}

echo $(date)": Ending jobs"
for node in $(scontrol show hostnames ${{SLURM_NODELIST}}); do
    srun -N 1 -n 1 -c 1 --overlap ${{WORKER_SCRIPTS}}/stop_worker.sh &
    sleep 2
done
exit
"""
def get_seconds(str_time, percent_off=0.10):
    days = str_time.split('-')[0]
    str_time = str_time.split('-')[-1]
    try:
       days = int(days)
    except ValueError as e:
        days = 0

    total_seconds = days*(60*60*24)
    str_time = str_time.split(":")
    convert = {0:1,1:60,2:60*60}
    for i, t in enumerate(str_time[::-1]):
        total_seconds += int(t)*convert[i]

    return round(total_seconds*(1-percent_off/100))
    

def submit_new_job(form):

    form['time_sec'] = get_seconds(form['time'])

    filename = f"new_nodes_{uuid.uuid1()}.sh"
    sbatch_cmd = f"sbatch --parsable {filename}"
    with open(filename, 'w') as script:
        script.write(job_template.format(**form))
    
    _stdout, _stderr, error = run_sh_command(sbatch_cmd, show_stdout=False)
    if error != 0:
        print(f"ERROR: failed to execute sbatch command: {sbatch_cmd}, {_stderr}, {_stdout}")
        return None
    else:
        form['job'] = _stdout
        return form


if __name__ == '__main__':
    condor_job_queue = get_condor_job_queue()
    print(condor_job_queue,"\n\n")
    slurm_workers = get_current_slurm_workers()
    print(slurm_workers,"\n\n")
    workers_needed = {'regular' : need_new_nodes(condor_job_queue, slurm_workers, 'regular'), 
                        'exvivo' : need_new_nodes(condor_job_queue, slurm_workers, 'exvivo')}
    print(workers_needed)

    for machine, needed in workers_needed.items():
        if needed >= 1:
            form = default_form[machine]
            form['num_nodes'] = needed
            form['machine'] = machine
            # Don't submit jobs for now, but it is working
            # ret = submit_new_job(form=form)


