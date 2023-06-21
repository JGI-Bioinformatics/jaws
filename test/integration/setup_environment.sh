#!/bin/bash -l
set -euxo pipefail

echo "Loading functions"
source "./test/integration/utils.sh"

echo "Loading default config values"
source "./test/integration/configs/default.sh"

echo "Loading deployment-specific config"
validate_vars "JAWS_DEPLOYMENT_NAME"
export JAWS_DEPLOYMENT_NAME=`echo $JAWS_DEPLOYMENT_NAME | awk '{print tolower($0)}'`
source "./test/integration/configs/deployments/$JAWS_DEPLOYMENT_NAME.sh"

echo "Loading site-specific config"
validate_vars "JAWS_SITE_NAME"
export JAWS_SITE_NAME=`echo $JAWS_SITE_NAME | awk '{print tolower($0)}'`
source "./test/integration/configs/sites/$JAWS_SITE_NAME.sh"

echo "Port setup"
# JAWS_DEPLOYMENT_NAME should be upper case for setting *_PORT per deployment
export JAWS_DEPLOYMENT_NAME=`echo $JAWS_DEPLOYMENT_NAME | tr '[:lower:]' '[:upper:]'`
eval export JAWS_SUPERVISOR_PORT=\$JAWS_SUPERVISOR_PORT_$JAWS_DEPLOYMENT_NAME
eval export JAWS_AUTH_PORT=\$JAWS_AUTH_PORT_$JAWS_DEPLOYMENT_NAME
eval export JAWS_REST_PORT=\$JAWS_REST_PORT_$JAWS_DEPLOYMENT_NAME
eval export JAWS_CROMWELL_PORT=\$JAWS_CROMWELL_PORT_$JAWS_DEPLOYMENT_NAME
export JAWS_DEPLOYMENT_NAME=`echo $JAWS_DEPLOYMENT_NAME | tr '[:upper:]' '[:lower:]'`


echo "Validating input variables"
REQUIRED_VARS="
JAWS_AUTH_PORT
JAWS_CROMWELL_PORT
JAWS_DB_HOST
JAWS_DB_PORT
JAWS_DB_PW
JAWS_DEFAULT_CONTAINER
JAWS_DEPLOYMENT_NAME
JAWS_GLOBUS_EP
JAWS_GLOBUS_HOST_PATH
JAWS_GROUP
JAWS_INSTALL_BASEDIR
JAWS_LOAD_PYTHON
JAWS_LOG_LEVEL
JAWS_MAX_RAM_GB
JAWS_REF_DATA_DIR
JAWS_REST_PORT
JAWS_RMQ_HOST
JAWS_RMQ_PORT
JAWS_RMQ_PW
JAWS_SCRATCH_BASEDIR
JAWS_SITE_NAME
JAWS_SUPERVISOR_PORT
JAWS_USERS_GROUP
"
validate_vars "$REQUIRED_VARS"

echo "Defining compound paths"
source "./test/integration/templates/all.sh"

echo "Installing into ${JAWS_CONFIG_DIR}"
echo "Creating paths and setting permissions"
FOLDERS="
JAWS_BIN_DIR
JAWS_CONFIG_DIR
JAWS_INSTALL_DIR
JAWS_LOGS_DIR
JAWS_SCRATCH_DIR
JAWS_SUPERVISOR_DIR
"

setup_dirs "$FOLDERS" "$JAWS_GROUP" 750