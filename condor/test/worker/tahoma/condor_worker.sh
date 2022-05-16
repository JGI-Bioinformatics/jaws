#!/bin/bash
# Pagurus task performance monitoring wrapper
# `venv` should be installed by deploying JAWS SITE
# (or pagurus should be found from `TAHOMA_PERFORMANCE_METRICS_SCRIPT` in the .gitlab.yml
# `--path` should be set by `TAHOMA_PERFORMANCE_METRICS_DIR` in the .gitlab.yml
source /tahoma/mscjgi/jaws-install/jaws-prod/site/bin/activate
pagurus \
--move \
--user $USER \
--path /tahoma/mscjgi/jaws_jtm/monitoring-runs \
--outfile condor_$(hostname)_$SLURM_JOB_ID.csv &

# Start a condor worker
. /tahoma/mscjgi/condor/env_worker.sh
condor_master &
sleep infinity