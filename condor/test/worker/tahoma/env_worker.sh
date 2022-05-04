# The script should be sourced by /bin/sh or similar
mkdir -p /home/svc-jtm-user/condor/worker_log/${HOSTNAME%%.*}/log
mkdir -p /home/svc-jtm-user/condor/worker_log/${HOSTNAME%%.*}/spool
mkdir -p /home/svc-jtm-user/condor/worker_log/${HOSTNAME%%.*}/execute

CONDOR_CONFIG="/tahoma/mscjgi/condor/worker.config"
export CONDOR_CONFIG
PATH="/tahoma/mscjgi/condor/condor/bin:/tahoma/mscjgi/condor/condor/sbin:$PATH"
export PATH
if [ "X" != "X${PYTHONPATH-}" ]; then
  PYTHONPATH="/tahoma/mscjgi/condor/condor/lib/python:$PYTHONPATH"
else
  PYTHONPATH="/tahoma/mscjgi/condor/condor/lib/python"
fi
export PYTHONPATH
