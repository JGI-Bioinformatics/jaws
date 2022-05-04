mkdir -p $SCRATCH/condor/$(hostname)
mkdir -p $SCRATCH/condor/$(hostname)/log
mkdir -p $SCRATCH/condor/$(hostname)/execute
export PATH="/global/cfs/projectdirs/jaws/condor/condor/bin":"/global/cfs/projectdirs/jaws/condor/condor/sbin":$PATH
export CONDOR_CONFIG="/global/cfs/projectdirs/jaws/condor/condor_config.worker.cori20"
