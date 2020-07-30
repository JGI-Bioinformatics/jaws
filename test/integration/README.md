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

### Cori (server: cori20.nersc.gov)

    Service          | dev   | staging | prod
    -----------------+-------+---------+------
    central-auth     | 3001  | 3002    | 3003
    central-rest     | 5001  | 5002    | 5003
    cromwell         | 50101 | 50102   | 50103
    supervisord-jaws | 64101 | 64102   | 64103
    supervisord-jtm  | 64111 | 64112   | 64113
