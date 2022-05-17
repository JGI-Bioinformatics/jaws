#!/bin/bash
now=$(date)
if [[ $(curl -s localhost:9118/api/v1/label/__name__/values) ]]; then
    echo "$now: ok"
else
    echo "$now: Something's wrong with HTCondor exporter. Kill and restart it."
    kill $(ps -aef | grep CondorExporter.py | awk '{print $2}')
    sh /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/start_exporter.sh
fi
exit $?
[svc-jtm-user@twf1 CondorExporter]$ more start_exporter.sh
#!/bin/bash
export CONDOR_CONFIG=/tahoma/mscjgi/condor/central.config
source /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/venv/bin/activate
export PYTHONPATH=/tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter:$PYTHONPATH
python /tahoma/mscjgi/condor/monitor/htcondor-monitor/CondorExporter/exporter/CondorExporter.py &
