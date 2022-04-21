#!/bin/bash -l
source /global/cfs/cdirs/jaws/condor/env_central.sh
source /global/cfs/cdirs/jaws/condor/pool-manager/venv/bin/activate
condor_pool_add -c /global/cfs/cdirs/jaws/condor/pool-manager/test/condor/jaws_condor/jaws_condor.ini -s /global/cfs/cdirs/jaws/condor/pool-manager/test/condor/jaws_condor/
site_config/cori_config.ini -d
