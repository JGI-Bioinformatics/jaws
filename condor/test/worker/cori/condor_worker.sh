#!/bin/bash
# Pagurus task performance monitoring wrapper
# `venv` should be installed by deploying JAWS SITE
# (or pagurus should be found from `CORI_PERFORMANCE_METRICS_SCRIPT` in the .gitlab.yml
# `--path` should be set by `CORI_PERFORMANCE_METRICS_DIR` in the .gitlab.yml

source /global/cfs/cdirs/jaws/jaws-install/jaws-prod/site/bin/activate
pagurus \
--move \
--user $USER \
--path  /global/cscratch1/sd/jaws_jtm/monitoring-runs \
--outfile condor_$(hostname)_$SLURM_JOB_ID.csv &

. /global/cfs/cdirs/jaws/condor/env_worker.sh
condor_master &

sleep infinity
