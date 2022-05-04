# JAWS Condor Pool Manager

JAWS Condor Pool Manager is to manager compute nodes for the Condor pool


## Start Condor Central

Installation directories
- Cori: /global/cfs/cdirs/jaws/condor
- Tahoma: /tahoma/mscjgi/condor

How to start the Condor Central 
- Cori
```
$ cd /global/cfs/cdirs/jaws/condor
$ source env.sh
$ condor_master
```

- Tahoma
```
$ cd /tahoma/mscjgi/condor/jaws/condor
$ source env.sh
$ condor_master
```

If it starts without issue, those services will be running
```
ps -aef | grep condor
svc-jtm+ 108388      1  0 Apr27 ?        00:00:12 condor_master
svc-jtm+ 108389 108388  0 Apr27 ?        00:04:59 condor_procd -A /home/svc-jtm-user/condor/central_log/log/procd_pipe -L /home/svc-jtm-user/condor/central_log/log/ProcLog -R 1000000 -S 60 -C 275072
svc-jtm+ 108390 108388  0 Apr27 ?        00:00:59 condor_shared_port -f -p 9618
svc-jtm+ 108395 108388  0 Apr27 ?        00:02:54 condor_collector -f
svc-jtm+ 108402 108388  0 Apr27 ?        00:01:39 condor_negotiator -f
svc-jtm+ 108404 108388  0 Apr27 ?        00:01:03 condor_schedd -f
```

* Note: if not, please check the logs under 
- Cori: /global/cscratch1/sd/jaws_jtm/jaws-condor/central_log
- Tahoma: /home/svc-jtm-user/condor/central_log

* To create a pool manually
- Tahoma
```
$ cd /tahoma/mscjgi/condor/jaws/condor
$ sbatch condor_worker_medium.job

(Once a compute node is allocated)

$ condor_status
Name               OpSys      Arch   State     Activity LoadAv Mem     ActvtyTime

slot1@t61          LINUX      X86_64 Unclaimed Idle	 0.000 385344  0+14:24:48

               Total Owner Claimed Unclaimed Matched Preempting Backfill  Drain

  X86_64/LINUX     1     0	 0         1	   0          0        0      0

         Total     1     0	 0         1	   0          0        0      0
```

* Test a Condor submission
- Tahoma
```
$ cd /tahoma/mscjgi/condor/jaws/condor
$ source env.sh
$ cd /tahoma/mscjgi/condor/jaws/condor/test
$ condor_submit sleep.jdf
Submitting job(s).
1 job(s) submitted to cluster 596.
$ condor_q


-- Schedd: svc-jtm-user@twf1 : <130.20.235.49:9618?... @ 05/04/22 13:04:22
OWNER        BATCH_NAME    SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
svc-jtm-user ID: 596      5/4  13:04      _      1      _      1 596.0

Total for query: 1 jobs; 0 completed, 0 removed, 0 idle, 1 running, 0 held, 0 suspended
Total for all users: 1 jobs; 0 completed, 0 removed, 0 idle, 1 running, 0 held, 0 suspended

```

* Check a job status (https://htcondor.readthedocs.io/en/latest/users-manual/managing-a-job.html)
- Tahoma
```
$ condor_history 596 -af JobStatus
4
```


* Condor Job Status Code (https://pages.cs.wisc.edu/~adesmet/status.html)
```angular2html
JobStatus in job ClassAds

0	Unexpanded	U
1	Idle	I
2	Running	R
3	Removed	X
4	Completed	C
5	Held	H
6	Submission_err	E

```
## JAWS Condor Pool Manager Installation Steps

Please see "Installing all services" in the project's main `README.md` file.

Or please clone the `condor` directory to your local file system and create your local Python virtual environment and do `pip install`.

ex) Tahoma setup
```
[svc-jtm-user@twf1 condor]$ source ../venv/bin/activate
(venv) [svc-jtm-user@twf1 condor]$ pip install --editable .
(venv) [svc-jtm-user@twf1 condor]$ which condor_pool_add
/tahoma/mscjgi/condor/pool-manager/venv/bin/condor_pool_add
(venv) [svc-jtm-user@twf1 condor]$ which condor_pool_remove
/tahoma/mscjgi/condor/pool-manager/venv/bin/condor_pool_remo
```



### Configuration

Parameters and commands for SLURM and HTCondor should be specified
in the `jaws_condor.ini` file.

Per each site, please refer to the *.ini files in `site_config`.

### Usage Example

```
source /tahoma/mscjgi/condor/env.sh && condor_pool_add -c jaws_condor.ini -s site_config/tahoma_config.ini -l $SCRATCH/jaws-condor/logs/condor_pool_add.log -d
```


- condor_pool_add.sh - 
```
#!/bin/bash -l
jaws_condor_root=/tahoma/mscjgi/condor
source $jaws_condor_root/pool-manager/venv/bin/activate
source $jaws_condor_root/env_central.sh
condor_pool_add -c $jaws_condor_root/pool-manager/jaws_condor.ini -s $jaws_condor_root/pool-manager/condor/jaws_condor/site_configs/tahoma_config.ini -l /tahoma/mscjgi/scrat
ch/condor/condor_pool_add.log -d
```


- condor_pool_remove.sh - 
```
#!/bin/bash -l
jaws_condor_root=/tahoma/mscjgi/condor
source $jaws_condor_root/pool-manager/venv/bin/activate
source $jaws_condor_root/env_central.sh
condor_pool_remove -c $jaws_condor_root/pool-manager/jaws_condor.ini -s $jaws_condor_root/pool-manager/condor/jaws_condor/site_configs/tahoma_config.ini -l /tahoma/mscjgi/sc
ratch/condor/condor_pool_remove.log -d
```


ex) Cori Cronjob

```
*/3 * * * * /global/cfs/cdirs/jaws/condor/pool-manager/condor_pool_add.sh > /dev/null 2>&1
*/10 * * * * /global/cfs/cdirs/jaws/condor/pool-manager/condor_pool_remove.sh > /dev/null 2>&1
```

ex) Tahoma Cronjob

```
*/3 * * * * /tahoma/mscjgi/condor/pool-manager/condor_pool_add.sh > /dev/null 2>&1
*/10 * * * * /tahoma/mscjgi/condor/pool-manager/condor_pool_remove.sh > /dev/null 2>&1
```