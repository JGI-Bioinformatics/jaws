# Integration Testing/Deployment on Cori

This directory contains scripts, which are used to deploy JAWS to the Cori system. Those scripts
are driven by the Gitlab CI/CD system, as specified in the .gitlab-ci.yml file in the root of
this repository.

## Architecture

As JAWS needs to be portable, one of the core tenets of this effort is to not depend on operating
system specific or site specific features. JAWS as a whole needs only Python3 to run, the rest
of the dependencies can be bootstrapped from there. As we can not rely on process supervision by
the operating system, [supervisord](https://www.supervisord.org) was chosen to fulfill this role.

Those instances are always running and can be controlled using the supervisorctl command. Typical
actions are starting and stopping a service. Access is controlled by a unique key.

     supervisord.conf
           |
    +--------------+      spawns services
    | supervisord  | -----------------------> jaws-site
    +--------------+          |
           |                  |-------------> jaws-central
           | controls         |
           |                  |-------------> other services
    +--------------+
    |   user/ci    |
    +--------------+

Every system needs two instances of supervisord, for privilege seperation between services and
user workloads: one for JAWS and one for JTM/Cromwell.


## Common Commands

To see this in action see .gitlab-ci.yml .

Start the supervisors. Only necessary once, after startup of the machine hosting the services: 

    collabsu jaws
    /tmp/jaws-supervisord/bin/supervisord -c /tmp/jaws-supervisord/supervisord-jaws.conf 
    logout
    collabsu jaws_jtm
    /tmp/jaws-supervisord/bin/supervisord -c /tmp/jaws-supervisord/supervisord-jtm.conf

Check the status of JAWS services:

    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jaws.conf status
    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jtm.conf status

Start the JAWS services:

    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jaws.conf start
    /tmp/jaws-supervisord/bin/supervisorctl -c /tmp/jaws-supervisord/supervisord-jtm.conf start

## Starting the gitlab-runner on Cori20
After a maintenance, it is very likely that the runner will need to be restarted
in order to accomplish this use the following steps:

    collabsu jaws
    cd $CFS/m342/jaws_runner/usr/bin
    nohup ./gitlab-runner run &

You can then check the UI on gitlab to see if the runner is up and working.

To do this, on the sidebar go to `Settings > CI/CD > Runners` and check if
the green dot is next to the cori20 runner.  


## Starting the gitlab-runner on lrc-services
Starting the runner is a bit different on LRC. Currently, there is no way to
impersonate a service user and start the service. Instead this must be done
using a user systemd. 

There is already a script present that starts up the services located here:  

    /global/home/groups-sw/lr_jgicloud/dev/jtm-service-dev.sh
    /global/home/groups-sw/lr_jgicloud/prod/jtm-service-prod.sh

To run the service you will want to run the following command:  

    systemctl start jawsdev

The above wil start the dev version of the runner. 
