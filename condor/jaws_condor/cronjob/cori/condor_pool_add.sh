#!/bin/bash -l
jaws_condor_root=/global/cfs/cdirs/jaws/condor
source $jaws_condor_root/pool-manager/venv/bin/activate
source $jaws_condor_root/env_central.sh
condor_pool_add -c $jaws_condor_root/pool-manager/jaws_condor.ini -s $jaws_condor_root/pool-manager/condor/jaws_condor/site_configs/cori_config.ini -d