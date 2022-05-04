#!/bin/bash
. /global/cfs/cdirs/jaws/condor/env_worker_highmem_jgi_shared.sh

#if [ "z$SLURM_PROCID" = "z0" ] ; then
#	~/bin/post2slack "Started $SLURM_JOBID with $SLURM_NNODES nodes"
#fi

mkdir -p $SCRATCH/condor/$(hostname)
mkdir -p $SCRATCH/condor/$(hostname)/log
mkdir -p $SCRATCH/condor/$(hostname)/execute

condor_master &
sleep 10000
