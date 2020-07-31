# Integration Testing/Deployment on PNNL Cascade

[Cascade Documentation](https://www.emsl.pnnl.gov/MSC/UserGuide/compute_resources/cascade_overview.html)

## Configuration

Host: gwf1.emsl.pnl.gov

Service User (JAWS): svc-jtm-manager

Service User (JTM): svc-jtm-user

Slurm Account: mscjgi

Parallel Filesystem (scratch): /dtemp/mscjgi

## Common Commands

To see this in action see .gitlab-ci.yml .

Start the supervisors. Only necessary once, after startup of the machine hosting the services:

    sudo -u svc-jtm-manager /tmp/jaws-supervisord/bin/supervisord -c /tmp/jaws-supervisord/supervisord-jaws.conf
    sudo -u svc-jtm-user /tmp/jaws-supervisord/bin/supervisord -c /tmp/jaws-supervisord/supervisord-jtm.conf

Check the status of JAWS services:

    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jaws.conf status
    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jtm.conf status

Start the JAWS services:

    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jaws.conf start
    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jtm.conf start


