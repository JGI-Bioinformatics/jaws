#!/bin/bash
export CONDOR_CONFIG=/global/cfs/cdirs/jaws/condor/condor_central.config
source /global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter/venv/bin/activate
export PYTHONPATH=/global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
python /global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter/exporter/CondorExporter.py &
