# HTCondor Backend for JAWS/Cromwell

This allows to execute JAWS runs using HTCondor which is a specialized workload management system for compute-intensive jobs created by the Center for High Throughput Computing in the Department of Computer Sciences at the University of Wisconsin-Madison (UW-Madison). Each JAWS site runs Condor Central service and condor worker pool(s) can be added to the Central via SLURM batch system. 


* Note that HTCondor doesn't support any functions to create worker pools so JAWS also has a service called `JAWS Condor Pool Manage`that maintains Condor worker pool using SLURM.  



## Cromwell Configuration for HTCondor

To HTCondor as Cromwell backend for JAWS, the below condor commands need to be added to the JAWS Cormwell configuration file.
- submit-docker = condor_submit
- check-alive = condor_q
- kill = condor_rm
- job-id-regex = "(?sm).*cluster (\\d+)..*"


The details can be found at `https://cromwell.readthedocs.io/en/stable/backends/HTcondor/`
And the production configuration files can also be found from the JAWS code repository under `condor/jaws_condor/test/cromwell`.




## JAWS Condor Pool Manager

JAWS Condor Pool Manager is to manager compute nodes for the Condor pool. The manager can keep a `min` number of nodes can be constantly running for JAWS. Further, it controls the `max` number of nodes in a pool so that JAWS does not overtake the SLURM system.

The worker pool is categorized into 4 types, `small, medium, large, xlarge`. The current supported types and the amount of resources provided per each site is the below.

The configuration file for each site can be found:
- Cori: /global/cfs/cdirs/jaws/condor/pool-manager/condor/jaws_condor/site_configs/cori_config.ini
- Tahoma: /tahoma/mscjgi/condor/pool-manager/condor/jaws_condor/site_configs/tahoma_config.ini
- JGI: /global/home/groups-sw/lr_jgicloud/condor/pool-manager/condor/jaws_condor/site_configs/jgi_config.ini


- Memory space (GB)

|                          | small      | medium   | large   | xlarge   |
|--------------------------|------------|----------|---------|----------|
| Cori @ NERSC             | n/a        | ~118   | n/a     | 118~1450 |
| Tahoma @ EMSL            | n/a        | ~364   | n/a     | 364~1480 |
| JGI @ LBL IT | ~46 | 46~236 | 236~492 | n/a      |



- CPUs

|                          | small | medium | large | xlarge  |
|--------------------------|-------|--------|-------|---------|
| Cori @ NERSC             | n/a   | 32     | n/a    | 36      |
| Tahoma @ EMSL            | n/a   | 36     | n/a   | 36 |
| JGI @ LBL IT | 32  | 32 | 32    | n/a     |


- Max walltime (hr)

|                          | small | medium | large | xlarge |
|--------------------------|-------|--------|-------|--------|
| Cori @ NERSC             | n/a   | 72     | n/a   | 168    |
| Tahoma @ EMSL            | n/a   | 48     | n/a   | 48     |
| JGI @ LBL IT | 72  | 72     | 72    | n/a    |


### Start Condor Central

The HTCondor installation directories per each site are
- Cori: /global/cfs/cdirs/jaws/condor
- Tahoma: /tahoma/mscjgi/condor
- JGI: /global/home/groups-sw/lr_jgicloud/condor

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

** Note: if not, please check the logs under  
- Cori: /global/cscratch1/sd/jaws_jtm/jaws-condor/central_log
- Tahoma: /home/svc-jtm-user/condor/central_log
- JGI: /global/scratch/users/jaws/condor/central_log

** Note: Each worker node's log
- Cori: /global/cscratch1/sd/jaws_jtm/jaws-condor/worker_log
- Tahoma: /home/svc-jtm-user/condor/worker_log
- JGI: /global/scratch/users/jaws/condor/worker_log

** Note: Each slurm job log
- Cori: /global/cscratch1/sd/jaws_jtm/jaws-condor/slurm
- Tahoma: /home/svc-jtm-user/condor/slurm
- JGI: /global/scratch/users/jaws/condor/slurm

** To create a pool manually
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

** Test a Condor submission
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

** Check a job status (https://htcondor.readthedocs.io/en/latest/users-manual/managing-a-job.html)
- Tahoma
Checking a running job
```angular2html
$ condor_q (job_id)
```
To get a detailed job status
```angular2html
$ condor_q (job_id) -l
```

To check a completed job
```
$ condor_history 596 -af JobStatus
4
```


** Condor Job Status Code (https://pages.cs.wisc.edu/~adesmet/status.html)
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
### JAWS Condor Pool Manager Installation Steps

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



#### Configuration

Parameters and commands for SLURM and HTCondor should be specified
in the `jaws_condor.ini` file.

Per each site, please refer to the *.ini files in `site_config`.


#### Usage Example

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



#### Cronjob per site

To enable the pool manager service, create cronjobs per site. 

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




## Misc info and Q$A

### Basic runtime section for a user WDL and the default values

- cpu/memory example:

```angular2html

    runtime {
        cpu: 1
        memory: "5G"
    }
```

- The interpretation of "cpu" and "memory" in the runtime section of the wdl has changed slightly for condor.  For instance, if before you requested memory=100GB cpu=64 node=1 and nwpn=2 you would expect each worker (or wdl-task) to get 50GB of memory and 32 cpus cores. In htcondor, each wdl-task will get 100GB and 64 cpus.  All other parameter keys like "nwpn" will be ignored. 
- The cpu/memory requirements should be specified per task. If not, the default cpu (cpu=1) and memory (mem=2GB) will be used for the task


### Q&A
- Advantages?
  - more standardization
  - The advantage of the condor backend is that only memory and cpu tags are required (time is recommended) making wdls more standard (no more poolname, nwpn)
shorter queue wait times
  - Condor also has the ability to automatically use various sizes of workers on a node. If a task that was using an entire node completes, it may be filled by a number of smaller tasks. Multiple users tasks will be able to use a single node which should decrease queue wait times.
no time-out errors
  - WDL tasks will get re-submitted by condor if they were put onto a node that didn't have enough time for the task to comple. Although the run could be re-submitted multiple times, probabilistically, it should be rare.  For long running tasks (i.e. tasks that have a long "time" specified in the runtime section) they will be preferentially put on a long running node; we are still working on the solution for long running tasks.

- So there's no need to have the same runtime parameters for each task?
  
  Correct

- So is there a dedicated set of nodes for jaws now, or will each of those tasks have to wait in a slurm queue?

  We are maintaining a small pool (with 4 nodes currently) constantly. This pool will expand and shrink according to the volume of the workload automatically by requesting more resources from slurm. In general some of the work should start almost instantly from the available nodes and then more nodes will be requested if work is sitting in the htcondor queue.

- What happens if the time parameter is or is not supplied? Time is now time for that task, not the whole pipeline? 

  We don't use time yet, but we are considering to use it to help schedule jobs on slurm nodes that won't timeout before the estimated time of the job (specified in runtime of wdl). Time may also make your wdl more portable if slurm is used as a backend by another user.

- So time is recommended because you believe it might be coming back in the future and you don't want that upcoming change to break our WDLs?
  
  Time is recommended for portability because if you share your WDL with someone that uses the SLURM backend (provided by Cromwell team) then time is required.
  Currently, if Condor runs a task and it runs out of time on the node, then Condor will automatically retry.  It will not show up as a retry in the Cromwell logs.  Condor workers request the maximum available time from SLURM and are reused for many tasks.  We do release unused nodes if there isn’t work to be done, although we maintain a minimum number.

- slurm vs jtm vs htcondor
  - SLURM: requests a resource per task 
  - JTM: requests specified number of node(s) with cpu/memory/time requirements and starts specified number of workers on it/them. The resources on a node are shared evenly among the workers. A worker pool for a workflow is exclusive.
  - HTCondor: runs on resource pools composed of compute node(s). It partitions the resource in realtime for scheduling user jobs. The main `requirements` for HTCondor scheduling are number of CPUs, memory and disk spaces (per task). Any condor pool can be shared among multiple workflow runs.
 
- Any resource limit for user task? Like max number of nodes?
  The current min and max pools sizes are 4~30
 
- Run states update
  We are working on adding new function to report the details on the task status.



## Troubleshooting for a user run

### Debug steps example
 
1. Get cromwell id from `jaws status <jaws_run_id>`
2. Open `cromwell.out`
ex) From Tahoma
```
$ vi /tahoma/mscjgi/jaws-install/jaws-prod/logs/cromwell.out
```
3. Search `f82e24d9-315d-46d3-9b1a-1316496cf26f`
4. Go to the last found
```
Check the content of stderr for potential additional information: /tahoma/mscjgi/scratch/jaws-prod/cromwell-executions/nmdc_metag/f82e24d9-315d-46d3-9b1a-1316496cf26f/call-qc/rqc.jgi_rqcfilter/5e911a32-3e04-4c06-bca7-baaa0f15d67f/call-rqcfilter/shard-0/execution/stderr.
2022-05-05 17:22:09,712 cromwell-system-akka.dispatchers.engine-dispatcher-73 INFO  - WorkflowManagerActor Workflow f82e24d9-315d-46d3-9b1a-1316496cf26f failed (during ExecutingWorkflowState): Job jgi_rqcfilter.rqcfilter:0:1 exited with return code 1 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.
```
5. look at the stderr 

```
[svc-jtm-user@twf1 logs]$ less /tahoma/mscjgi/scratch/jaws-prod/cromwell-executions/nmdc_metag/f82e24d9-315d-46d3-9b1a-1316496cf26f/call-qc/rqc.jgi_rqcfilter/5e911a32-3e04-4c06-bca7-baaa0f15d67f/call-rqcfilter/shard-0/execution/stderr
```


7. Found java assertion error
```
Exception in thread "main" java.lang.AssertionError: #Class     Reads   Bases   ReadPct BasePct Notes
Input   123273338       18614274038     100.000 100.000
Output  121274880       18177407513     98.379  97.653
Duplicate       0       0       0.000   0.000
LowQuality      1122518 169239244       0.911   0.909
PolyG   204099  1082140 0.166   0.006   SubsetOfLowQuality
N       312202  46820425        0.253   0.252   SubsetOfLowQuality
Artifact        162436  24313708        0.132   0.131
Spike-in        0       0       0.000   0.000
ForceTrim       0       123273338       0.000   0.662
Adapter 680278  115105954       0.552   0.618
SipMap  0       0       0.000   0.000
ChloroMap       0       0       0.000   0.000
MitoMap 0       0       0.000   0.000
RiboMap 0       0       0.000   0.000
RiboKmer        0       0       0.000   0.000
Human   0       0       0.000   0.000
Mouse   0       0       0.000   0.000
Cat     0       0       0.000   0.000
Dog     0       0       0.000   0.000
Microbe 31318   4648287 0.025   0.025
Other   0       0       0.000   0.000


trr=1996550
ri=123273338
ri=123273338
ro=121274880
        at jgi.RQCFilterStats.toString(RQCFilterStats.java:98)
        at jgi.RQCFilterStats.toString(RQCFilterStats.java:92)
        at jgi.RQCFilter2.dehumanize(RQCFilter2.java:2658)
        at jgi.RQCFilter2.process(RQCFilter2.java:1389)
        at jgi.RQCFilter2.main(RQCFilter2.java:73)

```



### JAWS errors command

```angular2html
jaws errors 33919
```



### How to get the node name that is running a task

if you want to get into a compute node by `srun`

1. Get the condor job id for a task which is running
```
2022-05-05 08:47:57,812 cromwell-system-akka.dispatchers.backend-dispatcher-5406 
INFO  - DispatchedConfigAsyncJobExecutionActor [UUID(f82e24d9)nmdc_metag.stage:NA:1]: job id: 837
```
2. Get the node info by running 

```condor_q <condor_job_id> -af RemoteHost```


ex)
```angular2html
[svc-jtm-user@twf1 test]$ condor_submit sleep.jdf
[svc-jtm-user@twf1 test]$ condor_q


-- Schedd: svc-jtm-user@twf1 : <130.20.235.49:9618?... @ 05/06/22 10:38:45
OWNER        BATCH_NAME    SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
svc-jtm-user ID: 905      5/6  10:38      _      1      _      1 905.0

Total for query: 1 jobs; 0 completed, 0 removed, 0 idle, 1 running, 0 held, 0 suspended
Total for all users: 1 jobs; 0 completed, 0 removed, 0 idle, 1 running, 0 held, 0 suspended
[svc-jtm-user@twf1 test]$ condor_q 905 -af RemoteHost
slot1_1@t97

```

3. Use `squeue` command to get the SLURM job id for the machine name

4. `srun` into the node
```angular2html
srun --overlap --jobid 123456  --pty /bin/bash
```

