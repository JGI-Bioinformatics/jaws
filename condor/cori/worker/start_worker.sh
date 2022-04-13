#!/bin/sh

echo "Starting worker on $(hostname)"

export CONDOR_INSTALL=/global/common/software/m3792/htcondor
# export CONDOR_INSTALL=/global/cfs/projectdirs/jaws/condor/condor
export PATH=$CONDOR_INSTALL/bin:$CONDOR_INSTALL/sbin:$PATH

# Cleanup old logs from that node
rm -rf $SCRATCH/condor/$(hostname)

# Create new log directories
mkdir -p $SCRATCH/condor/$(hostname)/log
mkdir -p $SCRATCH/condor/$(hostname)/execute
mkdir -p $SCRATCH/condor/$(hostname)/spool

case "$HOSTNAME" in
nid* | exvivo*)
    export CONDOR_CONFIG=/global/cfs/projectdirs/jaws/condor/condor_config.worker.cori20
    ;;
*)
    exit
    ;;

esac

# Start master with startd based on host
condor_master

sleep infinity
