#!/bin/bash
export CONDOR_CONFIG=/tahoma/mscjgi/condor/central.config
source /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/venv/bin/activate
export PYTHONPATH=/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
export LD_LIBRARY_PATH=/cluster/apps/python/Python-3.8.1/build/lib:$LD_LIBRARY_PATH
python /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/exporter/CondorExporter.py &
