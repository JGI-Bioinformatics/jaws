# Deployment
There are multiple scripts that handle the creation of directories, python
virtual environments, configuration files and supervisord wrapper scripts (shims)
for starting JAWS services. 

### Environment Variables and config directory
In order to deploy JAWS, you will need to set some environment variables. 
These variables are validated in the `setup_environment.sh` script. Within the
`config` directory there is a sites (eg jgi, perlmutter, dori) directory that
defines all the site specific variables and a `deployments` directory that
defines variables based on the deployment name (eg dev, staging, prod). The

#### defaults.sh

| Variable                         | Description                                          | defaults                                                                       |
|----------------------------------|------------------------------------------------------|--------------------------------------------------------------------------------|
| JAWS_DOCS_URL                    | url for the jaws documentation                       | https://jaws-docs.readthedocs.io/en/latest/                                    |
| JAWS_LOG_LEVEL                   | logger level (DEBUG, INFO, WARNING, ERROR, CRITICAL) | INFO                                                                           |
| JAWS_CENTRAL_HOST                | The url for the central server                       | http://jaws.lbl.gov"                                                           | 
| JAWS_RMQ_HOST                    | hostname for the rabbitmq server                     | rmq.lbl.gov                                                                    |
| JAWS_RMQ_PORT                    | port for the rabbitmq server                         | 5672                                                                           |
| JAWS_DB_HOST                     | hostname for the mysql database                      | mysql.lbl.gov                                                                  |
| JAWS_DB_PORT                     | port for the mysql database                          | 3306                                                                           |
| JAWS_CONTAINER_TYPE              | container type for workflows to use                  | docker                                                                         |
 | JAWS_DEFAULT_CONTAINER           | default container to pull to run workloads           | ubuntu@sha256:b5a61709a9a44284d88fb12e5c48db0409cfad5b69d4ff8224077c57302df9cf |
| JAWS_PERFORMANCE_METRICS_CLEANUP | clean after number of seconds                        | 600                                                                            |
| JAWS_LOAD_PYTHON                 | module command for loading python                    |                                                                                |
| JAWS_PYTHON                      | name of the python executable                        | python3                                                                        |
| JAWS_FAST_SCRATCH_DIR            | path of the fast scratch filesystem                  |                                                                                |
 | JAWS_BIG_SCRATCH_DIR             | path of the big scratch filesystem                   |
| JAWS_GITLAB_RUNNER               | path to the gitlab runner executable                 |                                                                                |
| JAWS_GITLAB_RUNNER_CONFIG        | path to the toml config file                         ||
| JAWS_RPC_SERVER_NUM_THREADS      | number of threads to run rpc server                  | 5                                                                              |
| JAWS_SITE_DNS_NAME               | DNS for the site                                     |                                                                                |
| JAWS_SUPERVISOR_NODAEMON         | run supervisor with nodaemon                         |                                                                                |
| JAWS_SETFACL                     | whether setfacl can be used on filesystem (1/0)      | 1                                                                              |
| JAWS_FILE_SYNC_DELAY_SEC         | delay in file sync                                   | 0                                                                              |


#### sites/{SITE}
The next directory contains several shell files that define the site specific
environment variables. The following are the sites that we deploy to: 
- assembly
- aws
- dori
- jgi
- nmdc
- perlmutter
- tahoma

Each site contains the same defined variables with a few exceptions. Here are the following
site specific variables: 

| Variable                        | Description                                               | defaults |
|---------------------------------|-----------------------------------------------------------|----------|
| JAWS_INSTALL_BASEDIR            | base directory where JAWS will be installed               |          |
| JAWS_GLOBUS_EP                  | UUID of the globus endpoint for that site                 |          |
| JAWS_GLOBUS_HOST_PATH           | path to the base directory of the globus endpoint         |          |
| JAWS_LOAD_PYTHON                | module command for loading python                         |          |
| JAWS_PYTHON                     | name of the python executable                             |          |
| JAWS_GROUP                      | group name that the jaws user is a part of                |          |
| JAWS_USERS_GROUP                | group name that all users of JGI are a part of            ||
 | JAWS_SCRATCH_BASEDIR            | path to the scratch directory                             |          |
| JAWS_REF_DATA_DIR               | location of the jaws reference database                   ||
| JAWS_MAX_RAM_GB                 | max ram in GB that a job can use at that site             ||
| JAWS_MAX_CPU                    | max cpu that a job can use at that site                   ||
| JAWS_PERFORMANCE_METRICS_SCRIPT | path to performance script                                ||
| JAWS_PERFORMANCE_METRICS_DIR    | path to the metrics directory                             ||
| JAWS_SUPERVISOR_PORT_PROD       | port for supervisor in production                         ||
| JAWS_AUTH_PORT_PROD             | port for JAWS Central auth in production                  ||
| JAWS_REST_PORT_PROD             | port for JAWS Central rest in production                  ||
| JAWS_CROMWELL_PORT_PROD         | port for Cromwell server in production                    ||
| JAWS_SUPERVISOR_PORT_STAGING    | port for supervisor in staging                            ||
| JAWS_AUTH_PORT_STAGING          | port for JAWS Central auth in staging                     ||
| JAWS_REST_PORT_STAGING          | port for JAWS Central rest in staging                     ||
| JAWS_CROMWELL_PORT_STAGING      | port for Cromwell server in staging                       ||
| JAWS_SUPERVISOR_PORT_DEV        | port for supervisor in dev                                ||
| JAWS_AUTH_PORT_DEV              | port for JAWS Central auth in dev                         ||
| JAWS_REST_PORT_DEV              | port for JAWS Central rest in dev                         ||
| JAWS_CROMWELL_PORT_DEV          | port for Cromwell server in dev                           ||
| JAWS_QUEUE_WAIT_LARGE           | sbatch command string to get queue information from Slurm ||


#### Templates
In the templates directory, there are several files that are used to generate the
configuration files for JAWS Site as well as writing the shims used by either supervisor,
systemd or other deployment tools that start the services. JAWS Site services can be started
using a variety of different methods. These include virtual environments, apptainer images, and
scrontab job submission. 

### deploy.sh
This is the main execution script that drives the deployment of JAWS. The pipeline steps for
deployment of JAWS Site generally follows these steps:

    1. Stop services
    2. Setup deployment environment
    3. Write configuration files, and shims. Set correct permissions.
    4. Start services

The deploy script will take two arguments, the launch method and the install method. By launc method,
we mean the method that is used to start/stop services. This can be either through supervisor, systemd or scrontab.
New methods can be added into the `deploy_tools` directory. Install method means how the services are installed. This
can be either through a virtual environment or just by pulling the most up to date container image and running
`apptainer instance start`. 

## Supervisor Architecture
As JAWS needs to be portable, one of the core tenets of this effort is to not depend on operating
system specific or site specific features. JAWS as a whole needs only Python3 to run, the rest
of the dependencies can be bootstrapped from there. As we can not rely on process supervision by
the operating system, [supervisord](https://www.supervisord.org) was chosen to fulfill this role.

Those instances are always running and can be controlled using the supervisorctl command. Typical
actions are starting and stopping a service. Access is controlled by a unique key.

     supervisord.conf
           |
    +--------------+      spawns services
    | supervisord  | -----------------------> rpc-server
    +--------------+          |
           |                  |-------------> run-daemon
           | controls         |-------------> transfer-daemon
           |                  |-------------> task-log
           |                  |-------------> etc...
    +--------------+
    |   user/ci    |
    +--------------+

Every system needs an instance of supervisord for each service.

### Supervisor Summary
Below we describe how supervisor related services are restarted for JAWS. 

The step by step instructions are in a [wiki document](https://code.jgi.doe.gov/advanced-analysis/jaws/-/wikis/Re-starting-JAWS-after-Maintenance) that outlines what to do to restart things after a scheduled maintenance of cori or jgi. Essentailly, there are re-start scripts that need to be run on CORI, JGI & Cascade for dev,staging & prod.    

The scripts can be seen under `jaws/test/integration/start_supervisor_services`

* start_central_services.sh
* start_jaws_cori_services.sh
* start_jaws_jgi_services.sh
* start_jtm_cori_services.sh
 

The steps you take after a scheduled maintenance or during a CI/CD deployment are roughly the same. During deployment, the supervisord is run in the .gitlab-ci.yml. You can go there to see these commands in action.

### Background of Restarting Services
Again, each instance (dev, staging, prod) will have its own supervisors. You will want to use the appropriate executable depending on which instance you want to work with.  

Starting the supervisors is only necessary once, after startup of the machine hosting the services. These are examples only; the paths to the executables can be seen in the re-start scripts listed above, or in the .gitlab-ci.yml.

    $ <command> <jaws_user>
    $ jaws-supervisord-<INSTANCE>/bin/supervisord -c jaws-supervisord-<INSTANCE>/supervisord-jaws.conf 

where the <command> is one of the following... 

LRC/DORI:
    - command: sudo -u <user> -i  
    - user: jaws

TAHOMA:
    - command: sudo -u <user> -i  
    - user: svc-jtm-manager 

PERLMUTTER:
    - sshproxy -c jaws
    - ssh -i ~/.ssh/jaws jaws@perlmutter-p1.nersc.gov


Note: For Tahoma you will need to ssh to the host `twf1.emsl.pnl.gov` before you attempt to change
into the user. If you are logged into tahoma, a simple `ssh twf1` should work fine. `svc-jtm-manager` maps to
`jaws` will be `svc-jtm-user`. 

### Check the status of JAWS services:

    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jaws.conf status
    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jtm.conf status

### Start the JAWS services:

    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jaws.conf start
    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jtm.conf start

Note: there exists two supervisord processes, one for jaws and one for jtm,  even if there are not two
separate jaws and jtm users in use at the deployment site.  

## Systemd

On the Dori cluster, we have access to systemd. When you run `deploy.sh` and specify `systemd` as the deploy tool,
our pipeline will create shim files to start the services as well as service files to enable them to run on systemd.
The pipeline will then take care of stopping and starting the services. We use systemd along with `apptainer instance`
to run these on Dori. 

## Scrontab
Unlike Cori, Perlmutter does not have any dedicated workflow nodes that one can ssh into and
start persistent services. Instead, persistent services need to be scheduled onto the system using Scrontab.  By
using the command `scrontab -e` you open a file that can be edited and include the configuration for your
long running job.

## Starting the gitlab-runner on Perlmutter
After a maintenance, it is very likely that the runner will need to be restarted. Perlmutter is a unique case where
there are no specialized workflow nodes and instead persistent services need to be started using scrontab which is a 
Slurm capability. In order to add services to Slurm you will want to edit the scrontab file (`scrontab -e`).

Perlmutter deployment includes the writing of a shim file for starting the gitlab runner
on Perlmutter.

For more info of scrontab, check out the NERSC documentation [here](https://docs.nersc.gov/jobs/workflow/scrontab/#scrontab-slurm-crontab)


## Starting the gitlab-runner on LRC

The LRC runner requires the python module to be loaded (`module load python`), before the runner is started.
This ensures that the correct python is present in the environment so that it does not run into missing library issues 
since it'll pick up the system python instead. 

    $ sudo -u jaws -i
    $ module load python
    $ cd /global/home/groups-sw/lr_jgicloud/jaws_ci_runner
    $ nohup ./usr/bin/gitlab-runner run -c configuration/config.toml &

## Starting the gitlab-runner on TAHOMA

Currently the gitlab runner is being managed by EMSL. To contact them, join the #emsl-jgi-coordination channel on the JGI slack. 


## File cleanup (cron)

The input files and cromwell-executions folders must be purged regularly by cron at each compute-site.

`inputs` folder:

    0 2 * * 0 find $SCRATCH/jaws-dev/inputs -mindepth 2 -mtime +14 -exec rm -rf {} \; 2>/dev/null
    0 2 * * 5 find $SCRATCH/jaws-staging/inputs -mindepth 2 -mtime +14 -exec rm -rf {} \; 2>/dev/null
    0 2 * * 6 find $SCRATCH/jaws-prod/inputs -mindepth 2 -mtime +14 -exec rm -rf {} \; 2>/dev/null

`cromwell-executions` (outputs) folder:

    0 5 * * 0 find $SCRATCH/jaws-dev/cromwell-executions -mindepth 2 -maxdepth 2 -type d -mtime +14 -exec rm -rf {} \; 2>/dev/null
    0 5 * * 5 find $SCRATCH/jaws-staging/cromwell-executions -mindepth 2 -maxdepth 2 -type d -mtime +14 -exec rm -rf {} \; 2>/dev/null
    0 5 * * 6 find $SCRATCH/jaws-prod/cromwell-executions -mindepth 2 -maxdepth 2 -type d -mtime +14 -exec rm -rf {} \; 2>/dev/null
