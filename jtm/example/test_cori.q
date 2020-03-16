#!/bin/bash -l
#SBATCH -c 1
#SBATCH --mem-per-cpu=1GB
#SBATCH -q genepool
#SBATCH -t 00:05:00
#SBATCH --job-name=jtm_static_worker
#SBATCH -A fungalp

module unload python
source ~/venv/bin/activate
for i in {1..1}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm-worker -l debug --jobid $SLURM_JOB_ID -wt static -t 00:05:00 -ct 0.250000 -tp _jtm_inner_request_queue.sulsj &
done
wait
