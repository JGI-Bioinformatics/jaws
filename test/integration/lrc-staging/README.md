# Integration Testing/Deployment on LBNL Lawrencium

This directory contains scripts, which are used to deploy JAWS to the LBNL Lawrencium system. Those scripts
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

    /tmp/jaws-supervisord-staging/bin/supervisord -c /tmp/jaws-supervisord-staging/supervisord-jaws.conf

Check the status of JAWS services:

    /tmp/jaws-supervisord-staging/bin/supervisorctl -c /tmp/jaws-supervisord-staging/supervisord-jaws.conf status

Start the JAWS services:

    /tmp/jaws-supervisord-staging/bin/supervisorctl -c /tmp/jaws-supervisord-staging/supervisord-jaws.conf start


## Starting the gitlab-runner on lrc-services

To run the service you will want to run the following command as the "jaws" user:

`/global/home/groups-sw/lr_jgicloud/jaws_ci_runner/usr/bin/gitlab-runner "run" "--config" "/global/home/groups-sw/lr_jgicloud/jaws_ci_runner/configuration/config.toml"``

There is only one runner for all deployments at this Site (staging/staging/dev).
