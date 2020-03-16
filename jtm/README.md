JGI Task Manager (JTM)
======================




- Note
   . Online documentation: https://jaws-docs.readthedocs.io/en/latest/index.html


JTM Manager
----

- It accepts user request and 

JTM workers with custom pool
---------
Create a manual worker
```
$ jtm-worker -tp <pool_name>
```

- Test example: this creates a manual JTM worker under "test" pool and a manager. "jtm-submit" sends a user command, "ls", to the worker to process.
```
$ jtm-manager &
$ jtm-worker -tp test
$ jtm-submit -cr "ls" -p test 
```

Create a set of workers

- If a task json has the field, "pool", JTM will create a pool of dynamic workers and send the tasks which have the same pool name.
- Dynamic workers will last only the specified job time.
- Example task json: the below requests a node from Cori HPC and starts 4 JTM workers on an allocated node. The workers will be run for 10 minutes and they will create a pool named as "test". 
```
{
    "command": "ps -aef | grep jtm > /tmp/jtm_ps_grep.out",
    "output_dir": "/tmp",
    "output_files": "/tmp/jtm_ps_grep.out",
    "pool": {"cluster": "cori",
             "poolname": "test",             
             "cpu": 32,
             "time": "00:10:00",
             "mem": "10GB",
             "node": 1,
             "nwpn": 4}
}
```

- Test example
```
$ jtm-manager 
$ jtm-submit -f example/jtm_task_pool.json
```

JTM User Interface
---------
- jtm-submit
- jtm-status 
- jtm-kill
- jtm-status
- jtm-isalive
- jtm-remove-pool
- jtm-check-manager
- jtm-check-worker
- jrm-resource-log

Job log and resource log
---------

- Raw resource consumption data file per task will be created under 'logs/resource'
- If a manager is started with "-r" flag, the manager's log also prints resource usage log from workers
ex)
>['child_pid', 'clone_time_rate', 'cpu_load', 'end_date', 'host_name', 'ip_address', 'job_time', 'life_left', 'mem_per_core', 'mem_per_node', 'num_cores', 'num_tasks', 'num_workers_on_this_node', 'perc_mem_used', 'pool_name', 'ret_msg', 'rmem_usage', 'root_pid', 'run_time', 'slurm_jobid', 'task_id', 'vmem_usage', 'worker_id', 'worker_type']
>83235,0.2,0.0,2018-11-15 14:24:21,ssul-dm.dhcp.lbl.gov,128.3.91.208,None,0,,,8,0,2,96.9,,hb,0.0,76199,0,0,637,0.0,kzrBimMXtPYxMk6ZXQrfWj,0
>83235,0.2,0.0,2018-11-15 14:24:28,ssul-dm.dhcp.lbl.gov,128.3.91.208,None,0,,,8,0,2,97.0,,hb,0.0,76199,0,0,637,0.0,kzrBimMXtPYxMk6ZXQrfWj,0


Test with Cromwell
=================

Start Cromwell server
-----
```
$ java -Dconfig.file=jtm.conf -jar cromwell-36.jar server
```

Start jtm
-----
```
$ jtm-manager -l debug
```

Start local jtm-worker
-----
```
$ jtm-worker -tp test -l debug
```

Start remote static jtm-workers
-----
- The flag, "-dr" is for dryrun. It will print the created batch script on screen.

```
$ jtm-worker -wt dynamic -cl cori -t 01:00:00 -m 10G -c 32 -N 1 -nwpn 1 -p test
```


