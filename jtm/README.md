JGI Task Manager (JTM)
======================



- Note
   . Online documentation: https://jaws-docs.readthedocs.io/en/latest/index.html


JTM Manager
----

- It accepts user task and processes the task on HPC

JTM workers with custom pool
---------
Create a manual worker (Environment variable 'JTM_CONFIG_FILE' can be used as well)
```
$ jtm --config=jtm.ini worker -p <pool_name>
```

- Test example: this creates a manual JTM worker under "test" pool and a manager. "jtm-submit" sends a user command, "ls", to the worker to process.
```
$ jtm --config=jtm.ini manager &
$ jtm --config=jtm.ini worker -p test
$ jtm --config=jtm.ini submit -cmd "ls" -p test 
```

Create a set of workers

- If a task json has the field, "pool", JTM will create a pool of dynamic workers and send the tasks which have the same pool name.
- Dynamic workers will last only the specified job time.
- Example task json: the below requests a node from Cori HPC and starts 4 JTM workers on an allocated node. The workers will be run for 10 minutes and they will create a pool named as "test". 
```
{
    "command": "ps -aef | grep jtm > /tmp/jtm_ps_grep.out",
    "pool": {"cluster": "cori",
             "poolname": "test",             
             "cpu": 32,
             "time": "00:10:00",
             "mem": "10GB",
             "node": 1,
             "nwpn": 1}
}
```


JTM User Interface
---------
```
Usage: jtm [OPTIONS] COMMAND [ARGS]...

Options:
  --debug
  --config TEXT  Config INI file
  --help         Show this message and exit.

Commands:
  check-manager  check if a JTM manager is running on a site :param ctx:...
  check-worker   total number of workers of the site if pool_name is...
  isalive        JtmInterface returns 0 if ready 1 if queued 2 if pending 3...
  kill           # JtmInterface returns code # 0: terminated successfully #...
  manager        JTM Manager Click wrapper :param ctx: :param log_dir:...
  remove-pool    Remove pool of workers from HPC This actually looks up...
  resource-log   Find resource log file from LOG_DIR for a task id and...
  status         JtmInterface returns 0 if ready 1 if queued 2 if pending 3...
  submit         JtmInterface returns 'task_id' if successfully queued jtm...
  worker         JTM Worker Click wrapper
```

Task resource usage log
---------

- Raw resource consumption data file per task will be created under 'logs/resource'
- If a manager is started with "-r" flag, the manager's log also prints resource usage log from workers
ex) 
```
$ jtm --config=jtm.ini manager -r
```

- Fields: "child_pid",  # 1
"clone_time_rate",  # 2
"cpu_load",  # 3
"end_date",  # 4
"host_name",  # 5
"ip_address",  # 6
"job_time",  # 7
"life_left",  # 8
"mem_per_core",  # 9
"mem_per_node",  # 10
"num_cores",  # 11
"num_tasks",  # 12
"num_workers_on_node",  # 13
"perc_mem_used",  # 14
"pool_name",  # 15
"ret_msg",  # 16
"rmem_usage",  # 17
"root_pid",  # 18
"run_time",  # 19
"slurm_jobid",  # 20
"task_id",  # 21
"vmem_usage",  # 22
"worker_id",  # 23
"worker_type",  # 24
"jtm_host_name",  # 25

ex)
```
RESOURCE : {1: 45992, 2: 0.0, 3: 88.3, 4: '2020-06-09 16:31:50', 5: 'sjsul-lm-2.local', 6: '127.0.0.1', 7: None, 25: 'ssul_laptop', 8: 0, 9: '', 10: '', 11: 8, 12: 0, 13: 1, 14: '99.0', 15: 'sulsj_test', 16: 'hb', 17: '432.6', 18: 45913, 19: 0, 20: 0, 21: 2839, 22: '122259.8', 23: 'wojPm5baWMNp7yRMGNgb8G', 24: 0, 26: 1}
RESOURCE : {1: 45992, 2: 0.0, 3: 97.3, 4: '2020-06-09 16:31:50', 5: 'sjsul-lm-2.local', 6: '127.0.0.1', 7: None, 25: 'ssul_laptop', 8: 0, 9: '', 10: '', 11: 8, 12: 0, 13: 1, 14: '99.7', 15: 'sulsj_test', 16: 'hb', 17: '700.4', 18: 45913, 19: 0, 20: 0, 21: 2839, 22: '164409.6', 23: 'wojPm5baWMNp7yRMGNgb8G', 24: 0, 26: 1}>83235,0.2,0.0,2018-11-15 14:24:28,ssul-dm.dhcp.lbl.gov,128.3.91.208,None,0,,,8,0,2,97.0,,hb,0.0,76199,0,0,637,0.0,kzrBimMXtPYxMk6ZXQrfWj,0
RESOURCE : {1: 45992, 2: 0.0, 3: 98.1, 4: '2020-06-09 16:31:50', 5: 'sjsul-lm-2.local', 6: '127.0.0.1', 7: None, 25: 'ssul_laptop', 8: 0, 9: '', 10: '', 11: 8, 12: 0, 13: 1, 14: '99.9', 15: 'sulsj_test', 16: 'hb', 17: '1098.0', 18: 45913, 19: 0, 20: 0, 21: 2839, 22: '206687.6', 23: 'wojPm5baWMNp7yRMGNgb8G', 24: 0, 26: 1}
```

Test with Cromwell
=================

Start Cromwell server
-----
```
$ java -Dconfig.file=jtm.conf -jar cromwell-VERSION.jar server
```

Start jtm
-----
```
$ jtm --config=jtm.ini manager 
```

Start local jtm-worker
-----
```
$ jtm --config=jtm.ini manager worker -p test
```

Start remote jtm-workers
-----
- The flag, "-dr" is for dryrun. It will print the created batch script on screen.

```
$ jtm --config=jtm.ini worker -wt dynamic -cl cori -t 01:00:00 -m 10G -c 32 -N 1 -nwpn 1 -p test
```

### Steps to install jtm
see the section "Installing all services" in the main README.md

