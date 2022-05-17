#!/bin/bash
export CONDOR_CONFIG=/global/home/groups-sw/lr_jgicloud/condor/condor_central.config
source /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter/venv/bin/activate
export PYTHONPATH=/global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
python /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter/exporter/CondorExporter.py &
