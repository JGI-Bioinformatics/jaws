# Deployment configs
## Overview
The configs directory contains a set of shell scripts and files to manage and configure various aspects of our
deployments across different environments and sites. Each file encapsulates a set of environment variables specific
to a deployment or a site, allowing for modularity and ease of management.

## Directory Structure
    .
    ├── README.md
    ├── default.sh
    ├── deployments
    │   ├── dev.sh
    │   ├── prod.sh
    │   └── staging.sh
    └── sites
        ├── assembly.sh
        ├── aws.sh
        ├── dori.sh
        ├── jgi.sh
        ├── nmdc.sh
        ├── perlmutter.sh
        └── tahoma.sh

### Key Components
#### default.sh
This shell script contains the default environment variables common to all deployments and sites. 
It usually includes settings that don't vary between environments, such as application-wide constants.

#### deployments
This directory contains shell scripts corresponding to different deployment environments. Each script sets variables specific to its deployment.

- **dev.sh**: Contains variables for the development environment.
- **prod.sh**: Contains variables for the production environment.
- **staging.sh**: Contains variables for the staging environment.

#### sites
This directory contains shell scripts for specific sites or platforms where the application is deployed. 
Each script sets variables unique to that particular site.
- **assembly.sh**: Site-specific settings for the Assembly environment.
- **aws.sh**: Site-specific settings for the AWS environment.
- **dori.sh**: Site-specific settings for the Dori environment.
- **jgi.sh**: Site-specific settings for the JGI environment.
- **nmdc.sh**: Site-specific settings for the NMDC environment.
- **perlmutter.sh**: Site-specific settings for the Perlmutter environment.
- **tahoma.sh**: Site-specific settings for the Tahoma environment.

### How to Use
1. Set `$JAWS_DEPLOYMENT_NAME` and `$JAWS_SITE_NAME` in your environment.
2. Setup Environment: Once those two variables are set, source `setup_environment.sh` which will set up all the environment variables for you
 
For example, if you are working in a development environment on the Tahoma site, you would source the variables in this order:

    export JAWS_DEPLOYMENT_NAME=dev
    export JAWS_SITE_NAME=tahoma
    source setup_environment.sh

## Environment Variables
In our setup, we have different shell scripts that are specific to the environment they are meant to operate in.
This document provides a detailed explanation of each environment variable used in the `dev.sh`, `staging.sh`,
and `prod.sh` scripts under the deployments directory.

### Environment Specific
#### dev.sh
`JAWS_LOG_LEVEL`  
Type: String  
Default: "DEBUG"  
Description: Sets the log level for the JAWS application when running in a development environment. Setting it to "DEBUG" will capture detailed log information useful for debugging.

`OCI_RUNTIME`
Type: String  
Default: Not set  
Description: Specifies the runtime for Open Container Initiative (OCI) compliant containers. Leaving it empty uses the system's default runtime.  

#### staging.sh
`JAWS_LOG_LEVEL`
Type: String  
Default: "DEBUG"  
Description: Configures the log level for the JAWS application in a staging environment. "DEBUG" is often used in staging to capture more detailed logs for testing.  

#### prod.sh
`JAWS_LOG_LEVEL`
Type: String  
Default: "INFO"  
Description: Defines the log level for the JAWS application in a production environment. Typically set to "INFO" to reduce the volume of log output, capturing only necessary information for operation monitoring.  

`JAWS_RPC_SERVER_NUM_THREADS`  
Type: Integer  
Default: 16  
Description: Sets the number of threads allocated for the JAWS RPC server. Increasing the thread count can improve performance by allowing the server to handle more requests concurrently.  

### Site Specific
#### Common Environment Variables

`JAWS_INSTALL_BASEDIR` Specifies the installation directory of JAWS.  

`JAWS_GLOBUS_EP` Globus Endpoint ID used for data transfer.  

`JAWS_GLOBUS_HOST_PATH` Host path where the Globus data will be stored or fetched.  

`JAWS_LOAD_PYTHON` The command to load the Python environment (often a module load command).  

`JAWS_PYTHON` The Python executable to be used.  

`JAWS_GROUP` Specifies the user group that owns the JAWS process.  

`JAWS_USERS_GROUP` Specifies the users group that should have access to JAWS.  

`JAWS_SCRATCH_BASEDIR` Specifies the base directory for temporary data storage.  

`JAWS_REF_DATA_DIR` The directory where the reference data is stored. 

`JAWS_MAX_RAM_GB` The maximum RAM in GB that are available on a node at a site.

`JAWS_MAX_CPU` The maximum number of CPU cores that are available on a node at a site.

`JAWS_SUPERVISOR_PORT_PROD`, `JAWS_AUTH_PORT_PROD`, `JAWS_REST_PORT_PROD`, `JAWS_CROMWELL_PORT_PROD`
Port numbers used by various services in the production environment.

`JAWS_SUPERVISOR_PORT_STAGING`, `JAWS_AUTH_PORT_STAGING`, `JAWS_REST_PORT_STAGING`, `JAWS_CROMWELL_PORT_STAGING`
Port numbers used by various services in the staging environment.

`JAWS_SUPERVISOR_PORT_DEV`, `JAWS_AUTH_PORT_DEV`, `JAWS_REST_PORT_DEV`, `JAWS_CROMWELL_PORT_DEV`
Port numbers used by various services in the development environment.

#### Site Specific

**assembly.sh**   
- `JAWS_SCRATCH_BASEDIR`: Different S3 paths are used based on deployment name (jaws-assembly-${JAWS_DEPLOYMENT_NAME} and jaws-site-${JAWS_DEPLOYMENT_NAME}).  

**dori.sh**  
- `JAWS_PERFORMANCE_METRICS_SCRIPT`: Specifies the script for performance metrics.
- `JAWS_PERFORMANCE_METRICS_DIR`: Directory to store performance metrics.

**jgi.sh**  

**nmdc.sh**
- `JAWS_ACCESS_GROUP`: Access group name for NMDC.
- `JAWS_PERLMUTTER`: Indicates that the script runs on the Perlmutter system.
- `TSIGPATH`: Path for the TSIG authentication key to update DNS.
- `JAWS_GITLAB_RUNNER` and `JAWS_GITLAB_RUNNER_CONFIG`: Configuration for the GitLab runner.

**perlmutter.sh**
- `JAWS_PERLMUTTER`: Indicates that the script runs on the Perlmutter system.
- `TSIGPATH`: Path for the TSIG authentication key to update DNS.

**tahoma.sh**
- `JAWS_BIG_SCRATCH_DIR`: Specifies a directory for large temporary data storage.
