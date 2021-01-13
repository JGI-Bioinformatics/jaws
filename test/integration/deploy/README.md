# Deployment

There are multiple scripts that handle the creation of directories, python
virtual environments, configuration files and supervisord wrapper scripts (shims)
for starting JAWS services.

### Environment Variables

In order to deploy JAWS, you will need to set some environment variables. This
is already taken care of in the `.gitlab-ci.yml` file. The following are
variables that must be set for the deployment to be successful:

* JAWS_CENTRAL_DB_PW: MySQL db password for JAWS central
* JAWS_CENTRAL_RMQ_PW: RabbitMQ password for JAWS central
* JAWS_SITE_CORI_DB_PW: MySQL db password for JAWS site on cori
* JAWS_SITE_CORI_RMQ_PW: RabbitMQ password for JAWS site on cori
* JAWS_CROMWELL_CORI_DB_PW: MySQL db password for cromwell on cori
* JAWS_JTM_CORI_DB_PW: MySQL db password for JTM on cori
* JAWS_SITE_LBNL_RMQ_PW: RabbitMQ password for JAWS site on JGI
* JAWS_CROMWELL_LBNL_DB_PW: MySQL db password for cromwell on JGI
* JAWS_JTM_LBNL_DB_PW: MySQL db password for JTM on JGI
* JAWS_JTM_RMQ_CORI_PW: RabbitMQ password for JTM on JGI
* JAWS_SUPERVISORD_PW: supervisor password


### Tool: deploy-jaws

This is the main execution script that drives the deployment of JAWS. The usage is as follows:

```bash
deploy-jaws [-d] <yaml_config> <service> [<service> ...] <output_install_dir>
```

Example:

```
deploy-jaws deploy-cori.yaml client site cromwell jtm /my/install/dir
```

The `deploy_jaws` bash script first creates a python environment for deployment within the `output_install_dir`, then executes a python script called `install-mgr` to install and deploy JAWS to the desired location.

Specifying the `-d` option allows one to install an instance of JAWS in `development mode` by using `pip install -e` option. Default uses `pip install`. The `-e` option enables editing of the installation code directly from the source code directory and having the changes take place in real time in the install directory.

The input `yaml_config` contains a yaml structure that defines the following for each component See `JAWS deployment configuration` for more information.

The `output_install_dir` defines where to install the JAWS deployment files. The location of the output config files and python virtual environments containing the JAWS code for a given component is configurable within the yaml template files.

If the `output_install_dir` does not exists, one will be created. If it exists, the following directories within the `output_install_dir` are deleted if found:

* configs:  jaws, jtm and cromwell config files
* deploy:  deployment scripts for installing components and creating shim files
* shims:  wrapper scripts used by supervisor to run jaws processes
* central: python virtual env containing jaws central code
* client: python virtual env containing jaws client code
* site: python virtual env containing jaws site code
* jtm: python virtual env containing jtm code
* supervisor (not deleted if found):  scripts and config files for running supervisor


### Tool: install-mgr

This python script handles the bulk of the installation and deployment of JAWS. The tool is written to be generic and configurable, using an input `yaml_config` file to define the components to install and how the installation and deployment is to be executed for each component. In summary, the yaml structure within `yaml_config` allows one to define config files for each component (if applicable), and deployment execution files for installing/deploying the component (if applicable).


### Deployment configuration file

The following yaml files are available for deploying an instance either on cori or jgi:

* deploy-cori.yaml: preconfigured to install one or more jaws services on cori. Default settings are set for the dev environment.
* deploy-jgi.yaml: preconfigured to install one or more jaws services on jgi. Default settings are set for the dev environment.
* deploy-cori-local.yaml: preconfigured to install one or more jaws services on cori as a local development instance. Default settings are set for the dev environment.

The structure of the yaml file is described as follows.

Below is an example of an input `yaml_config` used by `install-mgr`:

```yaml
services:
  - client:
      config_template: templates/client/jaws-client.config
      config_output: $JAWS_INSTALL_DIR/configs/jaws-client.conf
      deploy_script: templates/client/deploy-client
  - central:
      config_template: templates/central/jaws-central.config
      config_output: $JAWS_INSTALL_DIR/configs/jaws-central.conf
      deploy_script: templates/central/deploy-central
      environment:
        JAWS_DEPLOY_CENTRAL: 1

service_group:
  all: [client, central]

environment:
  JAWS_INSTALL_DIR:
  JAWS_SRC_DIR:
  JAWS_DEPLOY_NAME: dev
  JAWS_SITE_NAME: JGI
  JAWS_DEVELOPER_MODE: 0
  JAWS_IS_LOCAL_DEV_ENV: 0
```

Components to be installed are defined under the `services` section as a list. The values under `services` is a list available of jaws services to be installed. This allows services that are dependent on the installation of other services to be run in the correct order.

Each component should be given an unique name with the following keys:

* name of component (service).
* a config_template to use for creating the config file for that component (optional).
* a config_output for creating the config file for that component (optional; required if specifying config_template).
* a deployment script used for installing the component's code and generating the shim file for supervisor (optional).
* an optional environment section containing a dictionary of environment variables to defined for the given service.

The input service names to the `deploy-jaws` or `install-mgr` tool is defined under the `services` section of the `yaml` file. Additionally, a group name under the `service_group` section may be used as a shorthand for installing multiple services. Both service names a service group names may be combined as an input to the deployment tool.

Because the services are defined as a list, the deployment tool will execute the deployment scripts based on the order of this list for any `deploy_script` key is defined for that service.

The `services_group` section describes one or more names that can be used to define multiple service names as an input to the deployment tool.

The `environement` section describes one or more environment variables to define when running the deployment tool. The environment variables are used to create the jaws config files and running the jaws deployment scripts for each service.

Environment variables are defined under the 'environment' section at the root level or within the `services` section. If an external environment variable (set in a shell) is defined for the given name, the value of the external environment variable is used. Otherwise, the value defined in the yaml file is used. In addition, the value of each environment variable defined in the yaml file may contain other environment variable names either defined externally or within the yaml file itself. Substitutions will be done by the deployment tool.

The following environment variables are predefined when running `deploy-jaws` and cannot be changed with external environment variables:

```
  JAWS_INSTALL_DIR
  JAWS_SRC_DIR
  JAWS_DEVELOPER_MODE
  JAWS_IS_LOCAL_DEV_ENV
```


#### Specifying deployment environments

In order to specify a different deployment environment, such as `staging` or `prod`, specific environment variables need to be defined properly. The follow describes the environment variables needed for deploying `staging` or `prod`.

Example of deployment for staging:

```
export JAWS_DEPLOY_NAME=staging
export JAWS_CENTRAL_AUTH_PORT=3002
export JAWS_CENTRAL_REST_PORT=5002
export JAWS_CROMWELL_PORT=50102
export JAWS_SUPERVISOR_JAWS_PORT=64102
export JAWS_SUPERVISOR_JTM_PORT=64112
```

Example of deployment for prod:

```
export JAWS_DEPLOY_NAME=prod
export JAWS_CENTRAL_AUTH_PORT=3003
export JAWS_CENTRAL_REST_PORT=5003
export JAWS_CROMWELL_PORT=50103
export JAWS_SUPERVISOR_JAWS_PORT=64103
export JAWS_SUPERVISOR_JTM_PORT=64113
```

##### JAWS config files

The following JAWS config files are used, if the corresponding service is installed:

* jaws-central.conf (for installing central)
* jaws-client.conf (for installing client)
* jaws-jtm.conf (for installing site)
* jaws-site.conf  (for installing site)
* cromwell-shifter.conf ( (for installing site; cromwell config for using shifter)
* cromwell-singularity.conf ( (for installing site; cromwell config for using singularity)
* supervisord-config-jaws (for installing central and/or site)
* supervisord-config-jtm (for installing site)


#### JAWS deployment scripts

The JAWS deployment scripts for client, central and site perform the following operations:

* create python virtual environment
* pip install the code into the virtual environment
* create a shim file for executing the daemon process using `supervisor`

For the cromwell deployment script, the `cromwell.jar` is either copied from an existing location on the file system or retrieved using `wget`. The JAWS_CROMWELL_JAR_FILE environment variable are used to define how to get the jar files. If the value contains a url, then `wget` is used. Otherwise, the specified file location is used to copy the jar file to the installation directory. Similarly, the client deployment script either copies or uses wget to create the womtool.jar file in the install dir. The JAWS_CLIENT_WOMTOOL_JAR environment variable is used to define this.


##### Supervisor deployment script

The supervisor deployment script is configured to run after all the components are installed, configured and deployed. This script is responsible for starting/shutting down `supervisord` (where applicable) and starting/stopping processes using `supervisorctl` for both jaws and jtm processes (depending on which components are installed). Two separate supervisor config files are used, `supervisor-config-jaws` and `supervisor-config-jtm`.

The `supervisor-config-jaws` is setup to run the following JAWS processes, if installed:

* jaws-site-jtm-rpc
* jaws-site-central-rpc
* jaws-central-rest
* jaws-central-rpc
* jaws-central-auth

The `supervisor-config-jtm` is setup to tun the following JAWS processes, if installed:

* jaws-jtm
* jaws-cromwell

When the deployment script is executed, the following steps are performed:

* create python virtual environment and install `supervisor` using `pip install`.
* create bash scripts for running `supervisord` and `supervisorctl`
* if any of the JAWS processes defined in the `supervisor-config-jaws` are to be run, shutdown the `supervisord` process if found and start/restart `supervisord`.
* if `local-dev` deployment name is used and any of the processes defined in the `supervisor-config-jtm` are to be run, shutdown the `supervisord` process if found and start/restart `supervisord`. For any other deployment names used, it is assumed that the `supervisord` process is already running. The reasoning here is that the JTM processes need to be run as a different user than the `jaws` user, thus `supervisord` is not restarted.
* run `supervisorctl` to start all processes (where applicable).


#### Local deployment environment

To install a local JAWS instance on cori, run

```bash
deploy-jaws deploy-cori-local.yaml <one or more service to install> output_dir
```

In order to run jaws using the local instance, do the following

```bash
export JAWS_CLIENT_CONFIG=<output_dir>/configs/jaws-client.conf
source <output_dir>/client/bin/activate
```

Note that if more than one users are install a local instance on the same machine, the following environment variable values need to be unique per user:

```
  JAWS_CENTRAL_AUTH_PORT
  JAWS_CENTRAL_REST_PORT
  JAWS_CROMWELL_PORT
  JAWS_SUPERVISOR_JAWS_PORT
  JAWS_SUPERVISOR_JTM_PORT
```