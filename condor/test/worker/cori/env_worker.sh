mkdir -p /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log/${HOSTNAME%%.*}/log
mkdir -p /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log/${HOSTNAME%%.*}/spool
mkdir -p /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log/${HOSTNAME%%.*}/execute

CONDOR_CONFIG="/global/cfs/cdirs/jaws/condor/condor_worker.config"
export CONDOR_CONFIG
PATH="/global/cfs/cdirs/jaws/condor/condor/bin:/global/cfs/cdirs/jaws/condor/condor/sbin:$PATH"
export PATH
jaws_jtm@cori20:~/ssul/condor> more env_worker.sh
mkdir -p /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log/${HOSTNAME%%.*}/log
mkdir -p /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log/${HOSTNAME%%.*}/spool
mkdir -p /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log/${HOSTNAME%%.*}/execute

CONDOR_CONFIG="/global/cfs/cdirs/jaws/condor/condor_worker.config"
export CONDOR_CONFIG
PATH="/global/cfs/cdirs/jaws/condor/condor/bin:/global/cfs/cdirs/jaws/condor/condor/sbin:$PATH"
export PATH
