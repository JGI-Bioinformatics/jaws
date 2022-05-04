# The script should be sourced by /bin/sh or similar
mkdir -p /home/svc-jtm-user/condor/central_log

CONDOR_CONFIG="/tahoma/mscjgi/condor/central.config"
export CONDOR_CONFIG
PATH="/tahoma/mscjgi/condor/condor/bin:/tahoma/mscjgi/condor/condor/sbin:$PATH"
export PATH
if [ "X" != "X${PYTHONPATH-}" ]; then
  PYTHONPATH="/tahoma/mscjgi/condor/condor/lib/python:$PYTHONPATH"
else
  PYTHONPATH="/tahoma/mscjgi/condor/condor/lib/python"
fi
export PYTHONPATH