#### JAWS Monitoring

_ _ _

The JAWS monitoring process uses Prometheus to store metrics and Grafana web app to display the dashboard. There are two dashboards available. To access the dashboard, go to the <a href="https://grafana.jgi.lbl.gov">Grafana Dashboard</a> and select the appropriate dashboard.

The summarized  <a href="https://grafana.jgi.lbl.gov/d/S7d2rMPMk/jaws-status?orgId=1&refresh=1m">JAWS Status Dashboard</a> and the more comprehensive <a href="https://grafana.jgi.lbl.gov/d/u2LkNawGz/jaws-system-monitoring?orgId=1&refresh=1m&from=now-30m&to=now">JAWS System Monitor Dashboard</a>.

There are two components (python scripts) that are used to gather metrics for prometheus. The `jaws-prometheus` collector python script performs http and rabbitmq requests to each of the JAWS service to determine it's status. Additionally, a `site-monitor` python flask REST app reports metrics for a given site such as disk usage and supervisor status.

The `site-monitor` tool is installed on each JAWS site such as cori, lrc-services and jaws.lbl.gov, one per deployment environment (prod, staging). The `jaws-prometheus` tool is installed on jaws.lbl.gov, one per deployment environment (prod, staging).



#### Deployment

The following are instructions on how to restart the site-monitor daemon process for each site.

Restarting site monitor on Cori20.nersc.gov:

```
ssh cori20.nersc.gov
collabsu jaws
cd /global/cfs/projectdirs/jaws/jaws-monitor-[prod|staging]/supervisor

# Check status
bin/supervisorctl -c supervisord-monitor.conf status

# If a message like this appears, then it's likely that supervisord is not running.
http://localhost:xxxx refused connection

# To start supervisord if not running:
bin/supervisord -c supervisord-monitor.conf

# If supervisorctl reports that site monitor is not running, restart process:
bin/supervisorctl -c supervisord-monitor.conf restart jaws-[prod|staging]:site-monitor
```

Restarting site monitor on lrc-services.lbl.gov:

```
ssh user@lrc-services.lbl.gov
sudo -u jaws -i
cd /global/home/groups-sw/lr_jgicloud/jaws-install/jaws-monitor-[prod|staging]/supervisor

# Check status
bin/supervisorctl -c supervisord-monitor.conf status

# If a message like this appears, then it's likely that supervisord is not running.
http://localhost:xxxx refused connection

# To start supervisord if not running:
bin/supervisord -c supervisord-monitor.conf

# If supervisorctl reports that site monitor is not running, restart process:
bin/supervisorctl -c supervisord-monitor.conf restart jaws-[prod|staging]:site-monitor
```

Restarting site monitor on jaws.lbl.gov:

```
ssh user@lrc-services.lbl.gov
ssh jaws.lbl.gov
sudo -u jaws -i
cd /opt/jaws/jaws-monitor-[prod|staging]/supervisor

# Check status
bin/supervisorctl -c supervisord-monitor.conf status

# If a message like this appears, then it's likely that supervisord is not running.
http://localhost:xxxx refused connection

# To start supervisord if not running:
bin/supervisord -c supervisord-monitor.conf

# If supervisorctl reports that site monitor is not running, restart process:
bin/supervisorctl -c supervisord-monitor.conf restart jaws-[prod|staging]:site-monitor
```



The following are instructions on how to restart the jaws-prometheus daemon process on jaws.lbl.gov.

Restarting jaws prometheus collector:

```
ssh user@lrc-services.lbl.gov
ssh jaws.lbl.gov
sudo -u jaws -i
cd /opt/jaws/jaws-monitor-[prod|staging]/supervisor

# Check status
bin/supervisorctl -c supervisord-monitor.conf status

# If a message like this appears, then it's likely that supervisord is not running.
http://localhost:xxxx refused connection

# To start supervisord if not running:
bin/supervisord -c supervisord-monitor.conf

# If supervisorctl reports that site monitor is not running, restart process:
bin/supervisorctl -c supervisord-monitor.conf restart jaws-[prod|staging]:jaws-prometheus
```
