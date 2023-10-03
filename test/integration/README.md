# Deployment
## Overview
This document outlines the deployment strategy for our project, which involves a 
set of Bash scripts and configuration templates aimed to simplify the deployment process. 
The main driver script deploy.sh orchestrates the deployment by performing various tasks like writing configurations, 
setting up shims, and initiating services through specified launch tools (e.g., supervisor, systemd).

## Directory Structure
Here's an outline of the repository directory structure:
- README.md: General instructions and information about the deployment scripts.

### configs/
- **default.sh**: Default configuration script.
- **deployments**: Contains deployment-related configurations.
- **sites**: Site-specific configurations.

### deploy_tools/
Scripts related to deployment utilities:

- **scrontab.sh**: Script for deploying through cron.
- **supervisor.sh**: Script for deploying through Supervisor.
- **systemd.sh**: Script for deploying through systemd.

### templates/
Contains template files to write actual configurations during deployment:
Various .sh and .conf and .templ files to generate configurations.

## Utility Scripts
- **install_jaws_site.sh**: Installs the application site.
- **setup_environment.sh**: Sets up the necessary environment for the application. 
  These environment variables are in the configs/ directory. See README.md in there for more info.
- **utils.sh**: General utility functions used in the deployment.
- **write_configs.sh**: Writes configuration files based on the templates.
- **write_shims.sh**: Writes shim scripts for compatibility between different launch tools and install methods.

## Main Driver Script: deploy.sh
### Functionality
**Setting Flags**: It sets bash flags for strict error handling.  
**Source Utility Scripts**: It sources utility scripts that contain functions needed during deployment.  
**Deployment Function**: Contains a deploy function that does the following:  
Accepts arguments for the launch tool (e.g., supervisor, systemd) and the install method (e.g., venv, apptainer-sif, apptainer-docker).  
Calls functions to write configurations, write shims, and install the application.  

usage: 
```bash 
./deploy.sh [launch_tool] [install_method] [version]
```

#### Parameters
**launch_tool**: The launch tool you want to use for deployment (e.g., supervisor, systemd).  
**install_method**: The installation method to be used (e.g., venv, apptainer-sif, apptainer-docker).  
**version**: The version of the application you want to install.  

#### Example
If the launch tool is supervisor, and you want to spin down/up a service using 
an apptainer SIF image you would run: 

```bash
./deploy.sh supervisor apptainer-sif 1.2.3
```

#### Templates
In the templates directory, there are several files that are used to generate the
configuration files for JAWS Site as well as writing the shims used by either supervisor,
systemd or other deployment tools that start the services. JAWS Site services can be started
using a variety of different methods. These include virtual environments, apptainer images, and
scrontab job submissions. 

## Supervisor
### Overview
This section explains the Supervisor configurations for the JAWS project. 
Supervisor is used to control and monitor processes on UNIX-like operating systems. 
The configurations are defined in a set of .conf files.

### Architecture
As JAWS needs to be portable, one of the core tenets of this effort is to not depend on operating
system specific or site specific features. JAWS as a whole needs only Python3 to run, the rest
of the dependencies can be bootstrapped from there. As we can not rely on process supervision by
the operating system, [supervisord](https://www.supervisord.org) was chosen to fulfill this role.

Those instances are always running and can be controlled using the supervisorctl command. Typical
actions are starting and stopping a service. Access is controlled by a unique key.
                        
                      includes
    supervisord.conf -----------> supervisor.site.conf, supervisor.cromwell.conf, supervisor.pool-manager.conf, etc
           |
    +--------------+      spawns services (supervisorctl start)
    | supervisord  | ------------------------------------------> rpc-server
    +--------------+                             |
           |                                     |-------------> run-daemon
           | controls                            |-------------> transfer-daemon
           |                                     |-------------> task-log
           |                                     |-------------> etc...
    +--------------+
    |   user/ci    |
    +--------------+

Every system needs an instance of supervisord for each service.

### Supervisor configuration
There are four different configuration files that are required for starting up
services. The main `supervisor.conf` and three included configuration files
used to start up the "back-end" services of JAWS which include jaws-site, jaws-pool-manager,
Cromwell and HTCondor.

#### Included Configuration Files
supervisor.site.conf
Defines the following program processes:

- **runs**: Starts the run daemon.
- **transfers**: Starts the transfer daemon.
- **perf-metrics**: Starts the performance metrics daemon.
- **rpc-server**: Starts the RPC server.
- **task-log**: Starts the task log service.

All these programs are grouped under [group:jaws-site].

supervisor.cromwell.conf
Defines a program process:

- **cromwell**: Starts the Cromwell workflow engine.

This program is grouped under [group:jaws-cromwell].

supervisor.pool-manager.conf
Defines the following program processes:

- **pool-manager**: Manages pools of resources.
- **htcondor-server**: Starts the HTCondor server for job scheduling.
These programs are grouped under [group:jaws-pool-manager].


### Background of Restarting Services
Each instance (dev, staging, prod) will have its own supervisors. You will want to use the appropriate executable depending on which instance you want to work with.  

Starting the supervisors is only necessary once, after startup of the machine hosting the services. These are examples only; the paths to the executables can be seen in the re-start scripts listed above, or in the .gitlab-ci.yml.

    $ ${JAWS_INSTALL_BASEPATH}/jaws-install/<SITE>-[dev|staging|prod]/bin/supervisord -c ${JAWS_INSTALL_BASEPATH}/jaws-install/<SITE>-[dev|staging|prod]/supervisor.conf 

You will want to make sure you are logged in as the "jaws" user for that site. Not all sites have the service user
as jaws. For example, the jaws user on dori is `svc-jaws`, and the jaws user on tahoma is `svc-jtm-manager`.

Note: For Tahoma you will need to ssh to the host `twf1.emsl.pnl.gov` before you attempt to change
into the user. 

### Check the status of JAWS services:

    $ ${JAWS_INSTALL_BASEPATH}/jaws-install/<SITE>-[dev|staging|prod]/bin/supervisorctl status

### Start the JAWS services:

To start all the services (e.g. rpc-server, run-daemon, jaws-cromwell)

    $ ${JAWS_INSTALL_BASEPATH}/jaws-install/<SITE>-[dev|staging|prod]/bin/supervisorctl start all

### Important Notes
Make sure to restart Supervisor after making any changes to the configuration files.
Double-check permissions for accessing log and pidfile locations.

## Apptainer on Tahoma, LRC, and Dori
### Overview
Apptainer is a containerization solution we use to package and deploy our services across multiple platforms, 
including Tahoma, LRC, and Dori. This document outlines how we use Apptainer in these environments and 
how we work around certain limitations such as Supervisor's inability to handle background processes. 
We also discuss the usage of shims to facilitate these deployments.

### Apptainer Run on Tahoma and LRC
On Tahoma and LRC, we utilize apptainer run to initiate our services. These services are stateless, 
do not require background processes, and exit after performing their designated tasks. 
The apptainer run command ensures that the environment is isolated and that services are executed in a clean, 
controlled setting.

Example usage:

    $ apptainer run --cleanenv --env-file="$JAWS_CONFIG_DIR/site.env" --no-home --mount src="${JAWS_CONFIG_DIR}/jaws-site.conf",dst=/etc/config/site/jaws-site.conf,ro --mount src="${JAWS_LOGS_DIR}",dst=/var/log/site "${JAWS_BIN_DIR}/site-${JAWS_SITE_VERSION}.sif" --log "/var/log/site/${SERVICE}.log" --log-level "${JAWS_LOG_LEVEL}" "${SERVICE}"

### Apptainer Instance on Dori
On Dori, we use apptainer instance because Supervisor can't handle background processes, 
which are essential for certain services. An Apptainer instance allows us to run background tasks effectively. 
Unlike apptainer run, the instance command keeps the container running in the background, 
allowing long-running processes to operate over an extended period. Apptainer instance also provides more functionality
than run such a monitoring, logging and entering the shell of an already running instance (`apptainer shell`). For
more information, consult the Apptainer [documentation](https://apptainer.org/docs/user/main/running_services.html#instances-running-services).

### Shim Scripts
#### Apptainer Run Shim
The shim is used for initiation services using `apptainer run` on Tahoma and LRC:

    #!/bin/bash -l
    apptainer run --cleanenv --env-file="$JAWS_CONFIG_DIR/site.env" --no-home --mount src="${JAWS_CONFIG_DIR}/jaws-site.conf",dst=/etc/config/site/jaws-site.conf,ro --mount src="${JAWS_LOGS_DIR}",dst=/var/log/site "${JAWS_BIN_DIR}/site-${JAWS_SITE_VERSION}.sif" --log "/var/log/site/${SERVICE}.log" --log-level "${JAWS_LOG_LEVEL}" "${SERVICE}"

#### Apptainer Instance Shim

    #!/bin/bash
    apptainer instance start --cleanenv --env-file="$JAWS_CONFIG_DIR/site.env" --no-home --mount src="${JAWS_CONFIG_DIR}/jaws-site.conf",dst=/etc/config/site/jaws-site.conf,ro --mount src="${JAWS_LOGS_DIR}",dst=/var/log/site --pid-file ${JAWS_LOGS_DIR}/jaws-${SERVICE}-${JAWS_DEPLOYMENT_NAME}.pid "${JAWS_BIN_DIR}/site-${JAWS_SITE_VERSION}.sif" jaws-${SERVICE}-${JAWS_DEPLOYMENT_NAME} --log "/var/log/site/${SERVICE}.log" --log-level "${JAWS_LOG_LEVEL}" "${SERVICE}"

These shims take care of:
1. Loading environment variables from `.env` file.
2. Mounting required files and directories into the container.
3. Specifying the log files and log levels


## Site Specific Information
This section explains some of the difference between sites, specifically between Dori and Perlmutter since
those deployments differ from Tahoma and JGI. Tahoma and JGI are considered the "vanilla" installs that use
apptainer containers and supervisord to start up services.

## Dori
On the Dori cluster, we have access to systemd. When you run `deploy.sh` and specify `systemd` as the deploy tool,
our pipeline will create shim files to start the services as well as service files to enable them to run on systemd.
The pipeline will then take care of stopping and starting the services. We use systemd along with `apptainer instance`
to run these on Dori. 

### Running systemd
To run systemd, you will want to login as the `svc-jaws` user. If you want to start
the services, you will run: 

```
$ systemctl --user enable [SERVICE FILE]
$ systemctl --user start [SERVICE FILE]
```

service files are located in the `$HOME/config/systemd/` directory and are generated
by the CI/CD pipeline. 

If you are unsure of what services are under systemd, you can run

```
$ systemctl --user status
```

You can even grep this output to see what is running under each environment.
Here is an example: 

```
svc-jaws@ln010.jgi:/clusterfs/jgi/groups/dsi/homes/svc-jaws$ systemctl --user status | grep staging
           ├─task-log-staging.service
           │ ├─2769418 Apptainer instance: svc-jaws [jaws-task-log-staging]
           ├─perf-metrics-daemon-staging.service
           │ ├─1897716 Apptainer instance: svc-jaws [jaws-perf-metrics-daemon-staging]
           ├─run-daemon-staging.service
           │ ├─2769151 Apptainer instance: svc-jaws [jaws-run-daemon-staging]
           ├─rpc-server-staging.service
           │ ├─2769007 Apptainer instance: svc-jaws [jaws-rpc-server-staging]
           ├─task-logger-receive-staging.service
           │ ├─1897786 Apptainer instance: svc-jaws [jaws-task-logger-receive-staging]
           ├─transfer-daemon-staging.service
           │ ├─2769532 Apptainer instance: svc-jaws [jaws-transfer-daemon-staging]
```

The service names follow this scheme:  `[SERVICE]-[ENV].service`. As an example, the `rpc-server`
service would have this service file name: `rpc-server-staging.service`.


#### Example
If you want to start up the `rpc-server` service then you will run the following:

```
$ systemctl --user enable $HOME/config/rpc-server-dev.service
$ systemctl --user start $HOME/config/rpc-server-dev.service
```

## Perlmutter
This section outlines the process for deploying the JAWS project on the Perlmutter cluster. The deployment involves a few key components:

1. A scronjob script (starter.sh) that runs every 5 minutes to submit a job.
2. The `jaws-site-cronjob` script that sets up the environment and runs the main application using Supervisor.

### Deployment steps
#### Pre-requisites
Ensure that you have access to the Perlmutter cluster and that all the relevant directories and files are available.

#### SLURM scronjob
Before diving into the cron job and job script, it's essential to understand scronjob in the 
SLURM context. SLURM's scronjob allows you to schedule tasks similar to the Unix cron, but they are managed by the 
SLURM job scheduler. This is beneficial in cluster environments where resources are shared.

#### Scronjob (starter.sh)
The scronjob is a shell script that performs the following actions:

Checks if there's already an instance of the job running in the queue. If so, it will exit.
Checks for any scheduled maintenance on the cluster and adjusts the job submission accordingly.
Submits the job using sbatch, setting up various SLURM parameters including job name, dependencies, time, and resources.

#### jaws-site-cronjob
The job script takes care of:

Setting the dynamic DNS.
Initiating Supervisor, which then starts all the services defined in the Supervisor config files.

#### How to Setup:
This is automatically submitted by the starter.sh scronjob.

### Key Points
The starter.sh script is a SLURM scronjob scheduled to run every 5 minutes.
If the script detects cluster maintenance, it will adjust the job submission parameters accordingly.
The job script sets up the environment and initiates Supervisor to manage the services.
Both scripts have log paths specified for monitoring and debugging. If you want to submit to slurm
directly instead of through the pipeline, you will want to run the `starter.sh` script.

## Log files
### Supervisor
Supervisor enables us to monitor and control processes. One of the helpful
features of Supervisor is its ability to capture standard output and standard
error streams and redirect them to log files. This document will guide you on how to 
look up logs managed by Supervisor. 

### Location 
All logs generated by Supervisor for JAWS are located in the following directory:
`$JAWS_INSTALL_DIR/$JAWS_DEPLOYMENT_NAME/logs/`

**Note:** `$JAWS_INSTALL_DIR` and `$JAWS_DEPLOYMENT_NAME` are environment variables that points to the
installation directory. As an example here is what the installation path directory looks like:
`/clusterfs/jgi/groups/dsi/homes/svc-jaws/jaws-install/dori-prod/log`

#### Tail Logs in Real Time
If you want to tail logs in real time you can run the following command:

`tail -f [LOGFILE_NAME]`

For example: 

`tail -f rpc-server.log`

### Systemd
Systemd, the init system and service manager for Linux distributions, utilizes 
its logging system called journald. This system captures Syslog messages, kernel logs, 
initial RAM disk and early boot messages, as well as messages from running systemd services. 
Here's a brief guide to accessing these logs:

#### Viewing logs with journalctl
1. View All Logs: To view all the logs collected by journald, simply run:
`journalctl --user -xe`. This will open the logs in edit mode. You will need the user flag to access logs
that belong only to the `svc-jaws` user. 

2. Filter by Service: If you're interested in logs from a specific systemd service, you can filter the output. 
   For example, to view logs for a service named rpc-server-dev.service, run:
   `journalctl -u rpc-server-dev.service`
 
3. Live Tailing: To view logs in real time run: `journalctl -f`

#### Persistence
By default, journald logs are stored in a volatile storage (like /run/log/journal) and are lost upon reboot

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

Currently, the gitlab runner is being managed by EMSL. To contact them, join the #emsl-jgi-coordination channel on the JGI slack. 


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
