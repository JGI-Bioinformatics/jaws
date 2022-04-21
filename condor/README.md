# JAWS Condor Pool Manager

JAWS Condor Pool Manager is to manager compute nodes for the Condor pool

## Installation Steps

Please see "Installing all services" in the project's main `README.md` file.

## Configuration

Parameters and commands for SLURM and HTCondor should be specified
in the `jaws_condor.ini` file.

Per each site, please refer to the *.ini files in `site_config`.

### Usage

```
source /global/u1/j/jaws_jtm/ssul/condor/env_central.sh && condor_pool_add -c jaws_condor.ini -s ./site_config/cori_config.ini -d
```


- condor_pool_add.sh - 
```
#!/bin/bash -l
source /global/cfs/cdirs/jaws/condor/env_central.sh
source /global/cfs/cdirs/jaws/condor/pool-manager/venv/bin/activate
condor_pool_add -c jaws_condor.ini -s ./site_config/cori_config.ini -d
```