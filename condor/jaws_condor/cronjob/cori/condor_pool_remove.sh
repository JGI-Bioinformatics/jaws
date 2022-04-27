#!/bin/bash -l
source /global/cfs/cdirs/jaws/condor/pool-manager/venv/bin/activate
source /global/cfs/cdirs/jaws/condor/env_central.sh
condor_pool_remove -c /global/cfs/cdirs/jaws/condor/pool-manager/condor/jaws_condor/jaws_condor.ini -s /global/cfs/cdirs/jaws/condor/pool-manager/condor/jaws_condor/site_confi
gs/cori_config.ini -d