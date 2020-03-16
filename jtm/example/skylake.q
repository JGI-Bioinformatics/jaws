#!/bin/bash -l
#SBATCH -c 32
#SBATCH --mem=480GB
#SBATCH -C skylake
#SBATCH -A fungalp
#SBATCH -q genepool_special
#SBATCH -t 00:10:00
#SBATCH --job-name=jtm_dynamic_worker
#SBATCH -o jtm_dynamic_worker_vn5NRQ7DSTHuWGvr4Vkc5f.out
#SBATCH -e jtm_dynamic_worker_vn5NRQ7DSTHuWGvr4Vkc5f.err


module unload python
source ~/venv/bin/activate
for i in {1..1}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm-worker --jobid $SLURM_JOB_ID -cl cori -wt dynamic -t 00:10:00 -ct 0.200000 -tp test -nw 1 -C skylake  -c 32 -m 480GB  -wi vn5NRQ7DSTHuWGvr4Vkc5f${i} -A fungalp -q genepool_special &
    sleep 1
done
wait
