# The script should be sourced by /bin/sh or similar
mkdir -p /global/scratch/users/jaws/condor/central_log/log
mkdir -p /global/scratch/users/jaws/condor/central_log/spool
mkdir -p /global/scratch/users/jaws/condor/central_log/excute

CONDOR_CONFIG="/global/home/groups-sw/lr_jgicloud/condor/condor_central.config"
export CONDOR_CONFIG
PATH="/global/home/groups-sw/lr_jgicloud/condor/condor/bin:/global/home/groups-sw/lr_jgicloud/condor/condor/sbin:$PATH"
export PATH
if [ "X" != "X${PYTHONPATH-}" ]; then
  PYTHONPATH="/global/home/groups-sw/lr_jgicloud/condor/condor/lib/python:$PYTHONPATH"
else
  PYTHONPATH="/global/home/groups-sw/lr_jgicloud/condor/condor/lib/python"
fi
export PYTHONPATH
