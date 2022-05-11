# The script should be sourced by /bin/sh or similar
mkdir -p /tmp/jaws-condor/log
mkdir -p /tmp/jaws-condor/spool
touch /tmp/jaws-condor/spool/job_queue.log
mkdir -p /tmp/jaws-condor/config

CONDOR_CONFIG="/global/home/groups-sw/lr_jgicloud/condor/condor_worker.config"
export CONDOR_CONFIG
PATH="/global/home/groups-sw/lr_jgicloud/condor/condor/bin:/global/home/groups-sw/lr_jgicloud/condor/condor/sbin:$PATH"
export PATH
if [ "X" != "X${PYTHONPATH-}" ]; then
  PYTHONPATH="/global/home/groups-sw/lr_jgicloud/condor/condor/lib/python:$PYTHONPATH"
else
  PYTHONPATH="/global/home/groups-sw/lr_jgicloud/condor/condor/lib/python"
fi
export PYTHONPATH
