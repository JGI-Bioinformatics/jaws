# HTCondor Prometheous exporter 


This is a REST server that returns various performance metrics from a HTCondor Central service to REST requests.
By default, it listens to the port 9118.
This is used by a Prometheous service to collect metrics from a HTCondor service.

* NOTE: 
In JAWS, each site runs this exporter for the site's HTCondor Central. Each site should call REST to this exporter to get a JSON data of performance metrics. 



## JGI

1. cd /global/home/groups-sw/lr_jgicloud/condor/monitor/
2. git clone https://github.com/niclabs/htcondor-monitor
3. module load python 
4. cd /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter 
5. . venv/bin/activate && pip install -r requirements.txt


## How to start exporter

### JGI
```
$ cd cd /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter 
$ export PYTHONPATH=/global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter;$PYTHONPATH
$ python ./exporter/CondorExporter.py
Exporter listening on localhost:9118

```

### TAHOMA

```

(venv) [svc-jtm-user@twf1 CondorExporter]$ realpath .
/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter
(venv) [svc-jtm-user@twf1 CondorExporter]$ export PYTHONPATH=/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter;$PYTHONPATH
(venv) [svc-jtm-user@twf1 CondorExporter]$ echo $PYTHONPATH
/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter
(venv) [svc-jtm-user@twf1 CondorExporter]$ python ./exporter/CondorExporter.py
Exporter listening on localhost:9118


```



## How to request metrics

Command

```angular2html
curl -s localhost:9118/api/v1/label/__name__/values
```

Output example

```

(venv) [1060] jaws@lrc-services.lbl.gov:~/ssul/condor/monitor/htcondor-monitor/CondorExporter> curl -s localhost:9118/api/v1/label/__name__/values
# HELP process_virtual_memory_bytes Virtual memory size in bytes.
# TYPE process_virtual_memory_bytes gauge
process_virtual_memory_bytes 191229952.0
# HELP process_resident_memory_bytes Resident memory size in bytes.
# TYPE process_resident_memory_bytes gauge
process_resident_memory_bytes 25804800.0
# HELP process_start_time_seconds Start time of the process since unix epoch in seconds.
# TYPE process_start_time_seconds gauge
process_start_time_seconds 1652339788.81
# HELP process_cpu_seconds_total Total user and system CPU time spent in seconds.
# TYPE process_cpu_seconds_total counter
process_cpu_seconds_total 0.28
# HELP process_open_fds Number of open file descriptors.
# TYPE process_open_fds gauge
process_open_fds 7.0
# HELP process_max_fds Maximum number of open file descriptors.
# TYPE process_max_fds gauge
process_max_fds 262144.0
# HELP python_info Python platform information
# TYPE python_info gauge
python_info{implementation="CPython",major="3",minor="6",patchlevel="8",version="3.6.8"} 1.0
# HELP condor_slot_activity_idle Is this slot idle
# TYPE condor_slot_activity_idle gauge
condor_slot_activity_idle{address="10.0.45.0",machine="n0000.jgi0",slot="1"} 1.0
condor_slot_activity_idle{address="10.0.45.2",machine="n0002.jgi0",slot="1"} 1.0
condor_slot_activity_idle{address="10.0.45.1",machine="n0001.jgi0",slot="1"} 1.0
condor_slot_activity_idle{address="10.0.45.3",machine="n0003.jgi0",slot="1"} 1.0
# HELP condor_slot_activity_busy Is this slot busy
# TYPE condor_slot_activity_busy gauge
condor_slot_activity_busy{address="10.0.45.0",machine="n0000.jgi0",slot="1"} 0.0
condor_slot_activity_busy{address="10.0.45.2",machine="n0002.jgi0",slot="1"} 0.0
condor_slot_activity_busy{address="10.0.45.1",machine="n0001.jgi0",slot="1"} 0.0
condor_slot_activity_busy{address="10.0.45.3",machine="n0003.jgi0",slot="1"} 0.0
# HELP condor_slot_state_owner Is this slot in the owner state
# TYPE condor_slot_state_owner gauge
condor_slot_state_owner{address="10.0.45.0",machine="n0000.jgi0",slot="1"} 0.0
condor_slot_state_owner{address="10.0.45.2",machine="n0002.jgi0",slot="1"} 0.0
condor_slot_state_owner{address="10.0.45.1",machine="n0001.jgi0",slot="1"} 0.0
condor_slot_state_owner{address="10.0.45.3",machine="n0003.jgi0",slot="1"} 0.0
# HELP condor_slot_state_claimed Is this slot in the claimed state
# TYPE condor_slot_state_claimed gauge
condor_slot_state_claimed{address="10.0.45.0",machine="n0000.jgi0",slot="1"} 0.0
condor_slot_state_claimed{address="10.0.45.2",machine="n0002.jgi0",slot="1"} 0.0
condor_slot_state_claimed{address="10.0.45.1",machine="n0001.jgi0",slot="1"} 0.0
condor_slot_state_claimed{address="10.0.45.3",machine="n0003.jgi0",slot="1"} 0.0
# HELP condor_slot_state_unclaimed Is this slot in the unclaimed state
# TYPE condor_slot_state_unclaimed gauge
condor_slot_state_unclaimed{address="10.0.45.0",machine="n0000.jgi0",slot="1"} 1.0
condor_slot_state_unclaimed{address="10.0.45.2",machine="n0002.jgi0",slot="1"} 1.0
condor_slot_state_unclaimed{address="10.0.45.1",machine="n0001.jgi0",slot="1"} 1.0
condor_slot_state_unclaimed{address="10.0.45.3",machine="n0003.jgi0",slot="1"} 1.0
# HELP condor_job_state_idle Number of jobs on the idle state for a given cluster and submitter
# TYPE condor_job_state_idle gauge
# HELP condor_job_state_running Number of jobs on the running state for a given cluster and submitter
# TYPE condor_job_state_running gauge
# HELP condor_job_state_held Number of jobs on the held state for a given cluster and submitter
# TYPE condor_job_state_held gauge
# HELP condor_job_state_completed Number of jobs on the completed state for a given cluster and submitter
# TYPE condor_job_state_completed gauge
# HELP condor_job_avg_running_time_seconds Average running time for completed jobs for the specific cluster and submitter
# TYPE condor_job_avg_running_time_seconds gauge


```