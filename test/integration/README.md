# Integration Testing/Deployment

These directories contains scripts, which are used to deploy JAWS onto different sites. Those scripts
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

## Ports

    Service          | dev   | staging | prod
    -----------------+-------+---------+------
    central-auth     | 3001  | 3002    | 3003
    central-rest     | 5001  | 5002    | 5003
    cromwell         | 50101 | 50102   | 50103
    supervisord-jaws | 64101 | 64102   | 64103
    supervisord-jtm  | 64111 | 64112   | 64113


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

Note: there exists two supervisord processes, one for jaws and one for jtm,  even if there are not two
separate jaws and jtm users in use at the deployment site.


## Starting the gitlab-runner on Cori20

After a maintenance, it is very likely that the runner will need to be restarted
in order to accomplish this use the following steps:

    collabsu jaws
    cd $CFS/m342/jaws_runner/usr/bin
    nohup ./gitlab-runner run &

You can then check the UI on gitlab to see if the runner is up and working.

To do this, on the sidebar go to `Settings > CI/CD > Runners` and check if
the green dot is next to the cori20 runner.

## Starting the gitlab-runner on LRC

`/global/home/groups-sw/lr_jgicloud/jaws_ci_runner/usr/bin/gitlab-runner "run" "--config" "/global/home/groups-sw/lr_jgicloud/jaws_ci_runner/configuration/config.toml"``
