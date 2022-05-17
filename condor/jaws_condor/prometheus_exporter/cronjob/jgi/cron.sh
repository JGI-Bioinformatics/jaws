#!/bin/bash
now=$(date)
if [[ $(curl -s localhost:9118/api/v1/label/__name__/values) ]]; then
    echo "$now: ok"
else
    echo "$now: Something's wrong with HTCondor exporter. Kill and restart it."
    kill $(ps -aef | grep CondorExporter.py | awk '{print $2}')
    sh /global/home/groups-sw/lr_jgicloud/condor/monitor/htcondor-monitor/CondorExporter/start_exporter.sh
fi
exit $?
