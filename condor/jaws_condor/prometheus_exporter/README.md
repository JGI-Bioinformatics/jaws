# HTCondor Prometheus exporter 

This is a REST server that returns various performance metrics from a HTCondor Central service to REST requests.
This is used by a Prometheus service to collect metrics from a HTCondor service.

NOTE: 
- In JAWS, each site runs this exporter. Each JAWS Site should call a REST to this exporter to get performance metric data whenever JAWS Monitor requests them via AMQP RPC call.  
- By default, it listens to the port `9118`.
- The exporter code repo: https://github.com/niclabs/htcondor-monitor


## Per site setup

### JGI

```angular2html
1. cd /global/home/groups-sw/lr_jgicloud/condor/monitor/
2. git clone https://github.com/niclabs/htcondor-monitor
3. module load python 
4. cd /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter 
5. python -m venv venv
6. . venv/bin/activate && pip install -r requirements.txt
```
### TAHOMA

```angular2html
1. cd /tahoma/mscjgi/condor/monitor
2. git clone https://github.com/niclabs/htcondor-monitor
3. module load python 
4. cd /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter 
5. python -m venv venv
6. . venv/bin/activate && pip install -r requirements.txt
```

### CORI
```angular2html
1. cd /global/cfs/cdirs/jaws/condor/monitor
2. git clone https://github.com/niclabs/htcondor-monitor
3. module load python 
4. cd /global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter
5. python -m venv venv
6. . venv/bin/activate && pip install -r requirements.txt
```

## How to start exporter

### JGI
```
1. cd /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter 
2. . venv/bin/activate
3. export PYTHONPATH=/global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
4. export CONDOR_CONFIG=/global/home/groups-sw/lr_jgicloud/condor/condor_central.config
5. python ./exporter/CondorExporter.py &
```

### TAHOMA

```
1. cd /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter  
2. . venv/bin/activate
3. export PYTHONPATH=/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
4. export CONDOR_CONFIG=/tahoma/mscjgi/condor/central.config
5. python ./exporter/CondorExporter.py &
```

Ex)
```angular2html
[svc-jtm-user@twf1]$ cd /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter
[svc-jtm-user@twf1 CondorExporter]$ . venv/bin/activate
(venv) [svc-jtm-user@twf1 CondorExporter]$ realpath .
/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter
(venv) [svc-jtm-user@twf1 CondorExporter]$ export PYTHONPATH=/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
(venv) [svc-jtm-user@twf1 CondorExporter]$ echo $PYTHONPATH
/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter
(venv) [svc-jtm-user@twf1 CondorExporter]$ python ./exporter/CondorExporter.py
Exporter listening on localhost:9118
```


### CORI
```angular2html
1. cd /global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter
2. . venv/bin/activate
3. export PYTHONPATH=/global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
4. export CONDOR_CONFIG=/global/cfs/cdirs/jaws/condor/condor_central.config
5. python ./exporter/CondorExporter.py &
```

## How to request metrics

Command

```angular2html
curl -s localhost:9118/api/v1/label/__name__/values
```

Output example (JGI)

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


Output example (TAHOMA)
```angular2html
$ curl -s localhost:9118/api/v1/label/__name__/values
# HELP process_virtual_memory_bytes Virtual memory size in bytes.
# TYPE process_virtual_memory_bytes gauge
process_virtual_memory_bytes 395886592.0
# HELP process_resident_memory_bytes Resident memory size in bytes.
# TYPE process_resident_memory_bytes gauge
process_resident_memory_bytes 26644480.0
# HELP process_start_time_seconds Start time of the process since unix epoch in seconds.
# TYPE process_start_time_seconds gauge
process_start_time_seconds 1652720610.95
# HELP process_cpu_seconds_total Total user and system CPU time spent in seconds.
# TYPE process_cpu_seconds_total counter
process_cpu_seconds_total 1.17
# HELP process_open_fds Number of open file descriptors.
# TYPE process_open_fds gauge
process_open_fds 8.0
# HELP process_max_fds Maximum number of open file descriptors.
# TYPE process_max_fds gauge
process_max_fds 65536.0
# HELP python_info Python platform information
# TYPE python_info gauge
python_info{implementation="CPython",major="3",minor="8",patchlevel="1",version="3.8.1"} 1.0
# HELP condor_slot_activity_idle Is this slot idle
# TYPE condor_slot_activity_idle gauge
condor_slot_activity_idle{address="172.22.65.19",machine="t19",slot="1"} 1.0
condor_slot_activity_idle{address="172.22.65.20",machine="t20",slot="1"} 1.0
condor_slot_activity_idle{address="172.22.65.33",machine="t33",slot="1"} 0.0
condor_slot_activity_idle{address="172.22.65.34",machine="t34",slot="1"} 1.0
condor_slot_activity_idle{address="172.22.65.43",machine="t43",slot="1"} 1.0
# HELP condor_slot_activity_busy Is this slot busy
# TYPE condor_slot_activity_busy gauge
condor_slot_activity_busy{address="172.22.65.19",machine="t19",slot="1"} 0.0
condor_slot_activity_busy{address="172.22.65.20",machine="t20",slot="1"} 0.0
condor_slot_activity_busy{address="172.22.65.33",machine="t33",slot="1"} 0.0
condor_slot_activity_busy{address="172.22.65.34",machine="t34",slot="1"} 0.0
condor_slot_activity_busy{address="172.22.65.43",machine="t43",slot="1"} 0.0
# HELP condor_slot_state_owner Is this slot in the owner state
# TYPE condor_slot_state_owner gauge
condor_slot_state_owner{address="172.22.65.19",machine="t19",slot="1"} 0.0
condor_slot_state_owner{address="172.22.65.20",machine="t20",slot="1"} 0.0
condor_slot_state_owner{address="172.22.65.33",machine="t33",slot="1"} 0.0
condor_slot_state_owner{address="172.22.65.34",machine="t34",slot="1"} 0.0
condor_slot_state_owner{address="172.22.65.43",machine="t43",slot="1"} 0.0
# HELP condor_slot_state_claimed Is this slot in the claimed state
# TYPE condor_slot_state_claimed gauge
condor_slot_state_claimed{address="172.22.65.19",machine="t19",slot="1"} 0.0
condor_slot_state_claimed{address="172.22.65.20",machine="t20",slot="1"} 0.0
condor_slot_state_claimed{address="172.22.65.33",machine="t33",slot="1"} 0.0
condor_slot_state_claimed{address="172.22.65.34",machine="t34",slot="1"} 0.0
condor_slot_state_claimed{address="172.22.65.43",machine="t43",slot="1"} 0.0
# HELP condor_slot_state_unclaimed Is this slot in the unclaimed state
# TYPE condor_slot_state_unclaimed gauge
condor_slot_state_unclaimed{address="172.22.65.19",machine="t19",slot="1"} 1.0
condor_slot_state_unclaimed{address="172.22.65.20",machine="t20",slot="1"} 1.0
condor_slot_state_unclaimed{address="172.22.65.33",machine="t33",slot="1"} 0.0
condor_slot_state_unclaimed{address="172.22.65.34",machine="t34",slot="1"} 1.0
condor_slot_state_unclaimed{address="172.22.65.43",machine="t43",slot="1"} 1.0
# HELP condor_job_state_idle Number of jobs on the idle state for a given cluster and submitter
# TYPE condor_job_state_idle gauge
# HELP condor_job_state_running Number of jobs on the running state for a given cluster and submitter
# TYPE condor_job_state_running gauge
# HELP condor_job_state_held Number of jobs on the held state for a given cluster and submitter
# TYPE condor_job_state_held gauge
# HELP condor_job_state_completed Number of jobs on the completed state for a given cluster and submitter
# TYPE condor_job_state_completed gauge
# HELP condor_job_avg_running_time_seconds Average running time for completed jobs for the specific cluster an127.0.0.1 - - [16/May/2022 11:37:35] "GET /api/v1/label/__name__/values HTTP/1.1" 200 4157
d submitter
# TYPE condor_job_avg_running_time_seconds gauge

```