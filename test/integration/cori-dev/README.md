# Integration Testing/Deployment on Cori

## Common Commands

To see this in action see .gitlab-ci.yml .

Start the supervisors. Only necessary once, after startup of the machine hosting the services: 

    collabsu jaws
    /tmp/jaws-supervisord-dev/bin/supervisord -c /tmp/jaws-supervisord-dev/supervisord-jaws.conf 
    logout
    collabsu jaws_jtm
    /tmp/jaws-supervisord-dev/bin/supervisord -c /tmp/jaws-supervisord-dev/supervisord-jtm.conf

Check the status of JAWS services:

    /tmp/jaws-supervisord-dev/bin/supervisorctl -c /tmp/jaws-supervisord-dev/supervisord-jaws.conf status
    /tmp/jaws-supervisord-dev/bin/supervisorctl -c /tmp/jaws-supervisord-dev/supervisord-jtm.conf status

Start the JAWS services:

    /tmp/jaws-supervisord-dev/bin/supervisorctl -c /tmp/jaws-supervisord-dev/supervisord-jaws.conf start
    /tmp/jaws-supervisord-dev/bin/supervisorctl -c /tmp/jaws-supervisord-dev/supervisord-jtm.conf start

## Starting the gitlab-runner on Cori20
After a maintenance, it is very likely that the runner will need to be restarted
in order to accomplish this use the following steps:

    collabsu jaws
    cd $CFS/m342/jaws_runner/usr/bin
    nohup ./gitlab-runner run &

You can then check the UI on gitlab to see if the runner is up and working.

To do this, on the sidebar go to `Settings > CI/CD > Runners` and check if
the green dot is next to the cori20 runner.
