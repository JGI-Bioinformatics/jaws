# Deployment

There are multiple scripts that handle the creation of directories, python
virtual environments, configuration files and supervisord wrapper scripts (shims)
for starting JAWS services.

### Environment Variables

In order to deploy JAWS, you will need to set some environment variables. This
is already taken care of in the `.gitlab-ci.yml` file. The following are
variables that must be set for the deployment to be successful:

* JAWS_VERSION: version of JAWS
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



### deploy-jaws

This is the main execution script that drives the deployment of JAWS. The usage is as follows:

```bash
deploy-jaws [-d] <deploy_name> <json_profile_config> <output_install_dir>
```

The `deploy_jaws` script first creates a python environment for deployment within the `output_install_dir`, then executes a python script called `install-mgr` to install and deploy JAWS to the desired location.

Specifying the `-d` option allows one to install an instance of JAWS in `development mode` by using `pip install -e` option. This  enables editing of the installation code directly from the source code directory and having the changes take place in real time in the install directory.

The input `deploy_name` allows one to select the appropriate values for defining the entries in the JAWS config and deployment files. Although the `deploy_name` is arbitrary, for deployment we use `prod`, `staging` and `dev`. More information is described below under the `install-mgr` section.

The input `json_profile_config` contains a json structure that defines the following for each component:

* what component (service) to install.
* a config template to use for creating the config file for that component (optional).
* a deployment script used for installing the component's code and generating the shim file for supervisor (optional).
* a json template describing the values to substitute in the config and deployment script using jinja syntax (optional).

The list of available `json_profile_config` files is described under the _**JAWS deployment configuration**_ section. More detail regarding the structure of the json file is described in the _**install-mgr**_ section below.

The `output_install_dir` defines where to install the JAWS deployment files. The location of the output config files and python virtual environments containing the JAWS code for a given component is configurable within the json template files.

If the `output_install_dir` does not exists, one will be created. If it exists, the following directories within the `output_install_dir` are deleted if found:

* configs:  jaws, jtm and cromwell config files
* deploy:  deployment scripts for installing components and creating shim files
* shims:  wrapper scripts used by supervisor to run jaws processes
* central: python virtual env containing jaws central code
* client: python virtual env containing jaws client code
* site: python virtual env containing jaws site code
* jtm: python virtual env containing jtm code
* supervisor (not deleted if found):  scripts and config files for running supervisor



### install-mgr

This python script handles the bulk of the installation and deployment of JAWS. The tool is written to be generic and configurable, using an input `json_profile_config` file to define the components to install and how the installation and deployment is to be executed for each component. In summary, the json structure within `json_profile_config` allows one to define config files for each component (if applicable), and deployment execution files for installing/deploying the component (if applicable). The `install-mgr` tool supports jinja syntax  within the config and deployment execution files. See the [jinja](https://jinja.palletsprojects.com/en/2.11.x/) documentation for more details. The jinja variables are translated using a final deployment json file that contains the contents of the `json_profile_config` along with any `deploy_json_template` json files that are defined within the `component` section of the json structure. Additionally, each `deploy_json_template` may itself contain jinja syntax that will be translated using the json structure from the `json_profile_config` file.

Below is the structure of input `json_profile_config` used by `install-mgr`:

```json
{
    "install": {
        "meta": {
            "output_dir": "$INSTALL_DIR", (environment vars will be substituted)
            "src_dir": "$SRC_DIR",
            "deploy_name": "$DEPLOY_NAME",
            "developer_mode": "$DEVELOPER_MODE"
        },
        "components": {
            "some_name": {
                "add": 1, (a value of zero skips the component. larger values have lower execution order)
                      "deploy_json_template": "some_json_template_file" (optional)
                "deploy_script": "some_deploy_script", (optional)
                "config_template": "some_config_template_file", (optional)
                "output_config_name": "location of output config file", (optional)
            }
        }
    }
}

```

Environment variables may be defined anywhere within the json structure. The `install-mgr` tool will look for any string value that starts with a`$` and optionally a var name enclosed in `{}` brackets, like `${SOME_VAR}`. If found and the environment variable is defined outside of the tool, the json value will be substituted with the environment variable.

Components to be installed are defined under the `components` section. Each component should be given an arbritaray name with the following keys:

* add = 0|non-zero int. A zero indictates that the component should not be installed/processed. A non-zero value will trigger that component to be installed/processed. If dependencies exists between the installation of components, a greater value could be defined to the component that has a dependency to the previous one. The lower the add value, the later the execution order.

* config_template = name of config file with optional jinja syntax. Specify either a relative path to the `install-mgr` tool or a full path to the template file. This entry is optional.

* output_config_name = name of output config file. Any jinja syntax will be translated and the output is written to the file specified here. This entry is optional, though it is required if config_template is defined.

* deploy_script = name of script to execute for installing/deploying or configuring the given component. Optional jinja syntax may be defined in the file. This entry is optional.

* deploy_json_template = name of template json file defining the values for translating jinja variables within the config_template and/or deploy_script. This entry is optional. Note that any duplicate entries already defined in the
running deploy.json install file will be overwritten by the current entries in the json file.

Note that any additional key/value pairs defined in the json structure will be ignored by the `install-mgr` tool. However, additional json definitions in this file may be useful as they could be used to define the values within the `json_template_file`, `config_template_file` and/or `deployment_script` using jinja translation.

#### Specifying deployment environments

For installations that have multiple deployment environments (e.g., prod, staging, dev), the `install-mgr` supports a special json structure that allows one to define and use values for a given input deployment environment name. To define different values for each deployment name,  a `deploy_name` key must be defined with the sub keys corresponding to each deployment name within the `deploy_name` entry. Each sub key may have any json structure to define the values to use for jinja substitution. The `deploy_name` entry may be defined anywhere within the json structure inside either the `json_profile_config` or the `deploy_json_template` files. An example of the `deploy_name` structure is represented as follows.

```json
"site": {
    "deploy_name": {
        "dev": {
            "rabbitmq_port": 1,
            "mysql_port": 2
        },
        "prod": {
            "rabbitmq_port": 3,
            "mysql_port": 4
        }
    }
}
```

When running `install-mgr` with the input deployment name specified as `prod`, the json structure will be redefined as follows.

```json
"site": {
    "deploy_name": {
        "rabbitmq_port": 3,
        "mysql_port": 4
    }
}
```

An example of a config file that uses the `deploy_name` entry may look like this:

```json
rabbitmq_port={{site.deploy_name.rabbitmq_port}}
mysql_port={{site.deploy_name.mysql_port}}
```

#### Tool implementation

The following steps are performed by the `install-mgr` tool.

1. Translate any environment variables within the `json_profile_config` json file and write to output file in the installation directory.
2. For each component in which the `add` entry is non-zero, look for the `deploy_json_template` entry and if found, translate any jinja syntax within the defined file and merge the json structure to the running structure within the json file in the installation directory. Write updated json structure to the installation directory when done.
3. For each component in which the `add` entry is non-zero, look for the `config_template` entry and if found, translate any jinja syntax within the defined file and write the output to the name defined in the `output_config_name` entry. Jinja transation is performed using the final json file created in step 2.
4. For each component in which the `add` entry is non-zero, look for the `deploy_script` entry and if found, translate any jinja syntax within the defined file and write the output to the installation directory. Next, execute the script from within the installation directory.



### JAWS deployment configuration

For each site, a `json_profile_config` has been created to define which components are to be installed and where/how each component is installed and deployed.

* profile-cori.json: for deploying JAWS client, site, JTM and cromwell on cori.

* profile-jgi.json: for deploying  JAWS central on LabIT svm.

* profile-cori-local.json: for deploying all JAWS components on cori with a separate instance.

* profile-jgi-local.json: for deploying all JAWS components on jgi with a separate instance.

Here, a `component` corresponds to a JAWS daemon process. The available components are jaws-central-rest, jaws-central-auth, jaws-central-rpc, jaws-cromwell, jaws-jtm, jaws-site-central-rpc, jaws-site-daemon, jaws-site-jtm-rpc.

##### JAWS config files

The following JAWS config files are used, if component is installed:

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

For the cromwell deployment script, the `cromwell.jar` and `womtool.jar` is either copied from an existing location on the file system or retrieved using `wget`. The locations (file system or http) are defined in the `deploy_json_template` for cromwell using the `cromwell.jar_file` json entry and the `deploy_json_template` for client using the `client.womtool.jar` json entry, respectively. A shim file is created following the creation of the jar files.

##### Supervisor deployment script

The supervisor deployment script is configured to run after all the components are installed (`add=99`), configured and deployed. This script is responsible for starting/shutting down `supervisord` (where applicable) and starting/stopping processes using `supervisorctl` for both jaws and jtm processes (depending on which components are installed). Two separate supervisor config files are used, `supervisor-config-jaws` and `supervisor-config-jtm`.

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
deploy-jaws <deploy_name> profile-cori-local.json output_dir
```

where `deploy_name` may be `dev`, `staging`, or `prod`.

In order to run jaws using the local instance, do the following

```bash
export JAWS_CLIENT_CONFIG=<output_dir>/configs/jaws-client.conf
source <output_dir>/client/bin/activate
```

This installation will install all JAWS components (jaws-central-rest, jaws-central-auth, jaws-central-rpc, jaws-cromwell, jaws-jtm, jaws-site-central-rpc, jaws-site-daemon, jaws-site-jtm-rpc) within the specified `output_dir`. The following configuration describes the differences between a local deployment and a cori deployment environment aside from the type of components installed.

- central auth port: 3004
- central rest port: 5004
- client rest url: http://localhost
- client staging subdir `output_dir`/client/staging
- site cori jaws scratch_dir `output_dir`/scratch-`deploy_name`
- site cori globus uploads subdir: `output_dir`/uploads
- site cori globus down:w
loads subdir: `output_dir`/downloads
- jtm scratch dir: `output_dir`/scratch
- jtm worker config file: `output_dir`/config/jaws-jtm.conf
- jtm workker install dir: `output_dir`/jtm
- jtm environment activation: source `output_dir`/jtm-`deploy-name`/bin/activate
- jtm config file: `output_dir`/configs/jaws-jtm.conf
- cromwell docker root: `output_dir`/cromwell-executions
- cromwell worker log dir: `output_dir`/logs/cromwell-workflow-logs
- cromwell tmp dir: `output_dir`/cromwell-tmp
- cromwell port: 50104
- supervisor inet http server: 127.0.0.1
- supervisor url: http://localhost
- supervisor jaws port: 64104
- supervisor jtm port: 64114

This configuration is defined in the `templates/local-dev/cori-local.json` file.

