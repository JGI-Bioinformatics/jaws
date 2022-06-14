# Deployment

There are multiple scripts that handle the creation of directories, python
virtual environments, configuration files and supervisord wrapper scripts (shims)
for starting JAWS services. 

### Environment Variables 
In order to deploy JAWS, you will need to set some environment variables. This
is already taken care of in the `.gitlab-ci.yml` file. The following are
variables that must be set for the deployment to be successful:  

- DEPLOYMENT_NAME: Deployment environment name dev, staging, prod  
- JAWS_SITE: name of the deployment site ([SITE], JGI, CASCADE)  
- DEPLOY_JAWS_CLIENT: Whether client will be installed (1 or 0)  
- JAWS_SITES: All the different sites that JAWS is installed   
- JAWS_VERSION: version of JAWS  
- JAWS_DOCS_URL: URL where readthedocs is hosted ("https://jaws-docs.readthedocs.io/en/latest/")  
- JAWS_CENTRAL_HOST: hostname of JAWS central server  
- JAWS_RMQ_HOST: hostname of RabbitMQ server that all JAWS services will use
- JAWS_RMQ_PORT: port for RabbitMQ server that all JAWS services will use
- JAWS_GLOBUS_CLIENT_ID: Globus Confidential Application Client ID
- CROMWELL_VERSION: version of Cromwell
- JAWS_DB_HOST: hostname for the MYSQL database 
- JAWS_DB_PORT: port for the MYSQL database 

These are deployment specific. The deployment placeholder stands for either
dev, staging or prod.  
- [deployment]_LOG_LEVEL: Python logging level (DEBUG)  
- [deployment]_JAWS_SUPERVISOR_PORT: Supervisord port number  
- [deployment]_JTM_SUPERVISOR_PORT: Supervisord port number for JTM service  
- [deployment]_JAWS_AUTH_PORT: auth service port number  
- [deployment]_JAWS_REST_PORT: rest service port number  
- [deployment]_CROMWELL_PORT: cromwell port number  

The following changes per site and are prefixed with the site name:  
- [SITE]_INSTALL_BASEDIR: base directory where JAWS services will be installed  
- [SITE]_GLOBUS_EP: Globus endpoint id   
- [SITE]_GLOBUS_ROOT_DIR: Root path where globus will upload files   
- [SITE]_GLOBUS_DEFAULT_DIR: Default directory for globus upload  
- [SITE]_PYTHON: The name of the python executable (eg python, python3)  
- [SITE]_LOAD_PYTHON: The module command for loading python. Can be "" (ie - no modules)  
- [SITE]_JAWS_GROUP: Linux group name for JAWS  
- [SITE]_JTM_GROUP: Linux group name for JTM  
- [SITE]_CLIENT_GROUP: Linux group name for client. Group name is something all JGI users can use  
- [SITE]_JAWS_SCRATCH_BASEDIR: basename of SCRATCH directory for JAWS; "-{deployment_name}" will be appended
- [SITE]_JTM_SCRATCH_BASEDIR: basename of SCRATCH directory for JTM. Can be same as JAWS; "-{deployment_name}" will be appended
- [SITE]_JAWS_SW_BASEDIR: Base directory where JAWS code is located  
- [SITE]_JTM_SW_BASEDIR: Base directory where JTM code is located  
- [SITE]_REF_DATA_DIR: Location of pre-staged reference data that is mounted to each task's container
- [SITE]_FAST_SCRATCH_DIR: Path of fast storage local to worker nodes
- [SITE]_BIG_SCRATCH_DIR: Path of large storage local to worker nodes
- [SITE]_SUPERVISOR_HOST: Supervisor hostname for site  
- [SITE]_CONTAINER_TYPE: Container type the site uses (shifter, singularity etc)  
- [SITE]_CLUSTER_QOS: quality of service to use for a particular cluster  
- [SITE]_CLUSTER_PARTITION: Slurm partition name of the cluster  
- [SITE]_CLUSTER_ACCOUNT: Slurm cluster account to use for job submission  
- [SITE]_CLUSTER_CONSTRAINT: Cluster node type (eg - haswell, knl, skylake) 
- [SITE]_MAX_RAM_GB: Maximum ram for the cluster
- [SITE]_LOAD_JAVA: Module load command for Java  
- [SITE]_INPUTS_DIR: staging area for wdls, json files, zipped subworkflows files
- [SITE]_DOWNLOADS_DIR: output location of processed data sets. All under shared group permissions.


### deploy-jaws 
This is the main execution script that drives the deployment of JAWS. It first
sets up the environment by calling `define-envs` which looks for different
environment variables defined above. Before deployment it will call on  
`supervisorctl` to stop the services. Then it will perform a fresh install of
JAWS by generating the python virtual environments, configs and shims. It will
then restart the services with the new configuration and python installs. I


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


## Common Commands

### Summary
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

    $ logout

    $ <command> <jtm_user>
    $ jaws-supervisord-<INSTANCE>/bin/supervisord -c jaws-supervisord-<INSTANCE>/supervisord-jtm.conf

where the <command> is one of the following... 

CORI:
    - command: collabsu   
    - user: jaws | jaws_jtm

LRC:
    - command: sudo -u <user> -i  
    - user: jaws | jaws

TAHOMA:
    - command: sudo -u <user> -i  
    - user: svc-jtm-manager | svc-jtm-user


Note: For cascade you will need to ssh to the host `twf1.emsl.pnl.gov` before you attempt to change
into the user. If you are logged into tahoma, a simple `ssh twf1` should work fine. `svc-jtm-manager` maps to
`jaws` and `svc-jtm-user` maps to `jaws_jtm` so any services that need to run as the user that submits jobs
will be `svc-jtm-user`. 

### Check the status of JAWS services:

    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jaws.conf status
    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jtm.conf status

### Start the JAWS services:

    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jaws.conf start
    $ jaws-supervisord-<INSTANCE>/bin/supervisorctl -c jaws-supervisord-<INSTANCE>/supervisord-jtm.conf start

Note: there exists two supervisord processes, one for jaws and one for jtm,  even if there are not two
separate jaws and jtm users in use at the deployment site.  


## Starting the gitlab-runner on Cori20
(These instructions are also part of the re-start script described in [wiki document](https://code.jgi.doe.gov/advanced-analysis/jaws/-/wikis/Re-starting-JAWS-after-Maintenance).

After a maintenance, it is very likely that the runner will need to be restarted
in order to accomplish this use the following steps:

    $ collabsu jaws
    $ cd $CFS/m342/jaws_runner/usr/bin
    $ nohup ./gitlab-runner run &

You can then check the UI on gitlab to see if the runner is up and working.

To do this, on the sidebar go to `Settings > CI/CD > Runners` and check if
the green dot is next to the cori20 runner.

## Starting the gitlab-runner on LRC

The LRC runner requires the python module to be loaded (`module load python`), before the runner is started.
This ensures that the correct python is present in the environment so that it does not run into missing library issues 
since it'll pick up the system python instead. 

    $ sudo -u jaws -i
    $ module load python
    $ cd /global/home/groups-sw/lr_jgicloud/jaws_ci_runner
    $ nohup ./usr/bin/gitlab-runner run -c configuration/config.toml &

## Starting the gitlab-runner on Central
This is a special case since the gitlab runner on jaws.lbl.gov is a VM and is auto started by systemd as root.  Therefore, no intervention is required for re-starting it.
The below commands are what systemd runs:
    
    $ cd /usr/lib/gitlab-runner
    $ ./gitlab-runner run --working-directory /home/jaws --config /etc/gitlab-runner/config.toml --service gitlab-runner --syslog --user jaws

If necessary, use systemctl to check the status of the runner or to start/restart/stop the service:

    $ sudo systemctl [status|restart|stop|start] gitlab-runner.service

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
