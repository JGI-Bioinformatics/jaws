#!/bin/bash
if [[ $(curl -s localhost:9118/api/v1/label/__name__/values) ]]; then
    echo "ok"
else
    # Something's wrong with HTCondor exporter. Kill and restart it.
    kill $(ps -aef | grep CondorExporter.py | awk '{print $2}')
    sh /global/cfs/cdirs/jaws/condor/monitor/htcondor-monitor/CondorExporter/start_exporter.sh
fi
exit $?
