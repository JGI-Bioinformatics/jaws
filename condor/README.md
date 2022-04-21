# JAWS Condor Pool Manager

JAWS Condor Pool Manager is to manager compute nodes for the Condor pool

## Installation Steps

Please see "Installing all services" in the project's main `README.md` file.

Or please clone the `condor` directory to your local file system and create your local Python virtual environment and do `pip install`.
```

condor $> source ~/venv/bin/activate
(venv) condor $> pip install --editable .
...

(venv) condor $> which condor_pool_add
 ~/venv/bin/condor_pool_add

```

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


- condor_pool_remove.sh - 
```
#!/bin/bash -l
source /global/cfs/cdirs/jaws/condor/env_central.sh
source /global/cfs/cdirs/jaws/condor/pool-manager/venv/bin/activate
condor_pool_remove -c jaws_condor.ini -s ./site_config/cori_config.ini -d
```


ex) Cori Crontab

```
*/3 * * * * /global/cfs/cdirs/jaws/condor/pool-manager/test/condor/jaws_condor/condor_pool_add.sh >> /global/cscratch1/sd/jaws_jtm/condor/condor_pool_add.log 2>&1
*/10 * * * * /global/cfs/cdirs/jaws/condor/pool-manager/test/condor/jaws_condor/condor_pool_remove.sh >> /global/cscratch1/sd/jaws_jtm/condor/condor_pool_remove.log 2>&1
```