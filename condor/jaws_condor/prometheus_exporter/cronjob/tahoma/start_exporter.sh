#!/bin/bash
export CONDOR_CONFIG=/tahoma/mscjgi/condor/central.config
source /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/venv/bin/activate
export PYTHONPATH=/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
python /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/exporter/CondorExporter.py &
