#!/bin/bash -l

#set -eo pipefail

echo "BEGIN deploy-jaws-site"

echo "Loading functions"
source "./test/integration/utils.sh"

echo "Loading deployment-specific config"
validate_vars "JAWS_DEPLOYMENT_NAME"
export JAWS_DEPLOYMENT_NAME=`echo $JAWS_DEPLOYMENT_NAME | awk '{print tolower($0)}'`
source "./test/integration/configs/deployments/$JAWS_DEPLOYMENT_NAME.sh"

echo "Loading site-specific config"
validate_vars "JAWS_SITE_NAME"
export JAWS_SITE_NAME=`echo $JAWS_SITE_NAME | awk '{print tolower($0)}'`
source "./test/integration/configs/sites/$JAWS_SITE_NAME.sh"

echo "Validating input variables"
REQUIRED_VARS="
JAWS_AUTH_PORT
JAWS_CONTAINER_TMPDIR
JAWS_CONTAINER_TYPE
JAWS_CROMWELL_PORT
JAWS_DB_HOST
JAWS_DB_PORT
JAWS_DB_PW
JAWS_DEFAULT_CONTAINER
JAWS_DEPLOYMENT_NAME
JAWS_GLOBUS_CLIENT_ID
JAWS_GLOBUS_CLIENT_SECRET
JAWS_GLOBUS_EP
JAWS_GLOBUS_HOST_PATH
JAWS_GROUP
JAWS_INSTALL_BASEDIR
JAWS_LOAD_PYTHON
JAWS_LOG_LEVEL
JAWS_MAX_RAM_GB
JAWS_PERFORMANCE_METRICS_DIR
JAWS_PERFORMANCE_METRICS_SCRIPT
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
export JAWS_INSTALL_DIR="$JAWS_INSTALL_BASEDIR/$JAWS_DEPLOYMENT_NAME/$JAWS_SITE_NAME"
export JAWS_CONFIG_DIR="$JAWS_INSTALL_DIR/config"
export JAWS_BIN_DIR="$JAWS_INSTALL_DIR/bin"
export JAWS_LOGS_DIR="$JAWS_INSTALL_DIR/log"
export JAWS_SUPERVISOR_DIR="$JAWS_INSTALL_DIR/supervisor"
if [[ $JAWS_SCRATCH_BASEDIR =~ ^s3:\/\/ ]]; then
    # EACH AWS DEPLOYMENT HAS A SEPARATE S3 BUCKET
    export JAWS_SCRATCH_DIR="${JAWS_SCRATCH_BASEDIR}-${JAWS_DEPLOYMENT_NAME}"
else
    # EACH SERVER DEPLOYMENT HAS A SEPARATE NFS SUBDIR
    export JAWS_SCRATCH_DIR="$JAWS_SCRATCH_BASEDIR/jaws-$JAWS_DEPLOYMENT_NAME"
fi

# Folders for Run inputs and downloads from other jaws-sites
export JAWS_INPUTS_DIR="$JAWS_SCRATCH_DIR/inputs"
export JAWS_DOWNLOADS_DIR="$JAWS_SCRATCH_DIR/downloads"

# Cromwell-related dirs
export JAWS_CROMWELL_EXECUTIONS_DIR="$JAWS_SCRATCH_DIR/cromwell-executions"

# Performance metrics
export JAWS_PERFORMANCE_METRICS_SCRIPT="$JAWS_INSTALL_DIR/site/bin/pagurus"

# Globus has some rules about how paths are interpreted, depending on how the endpoint is configured.
[[ ${JAWS_GLOBUS_ROOT_DIR:-1} == "/" ]] || JAWS_GLOBUS_ROOT_DIR="$JAWS_GLOBUS_ROOT_DIR/"

# virtual envs
export JAWS_VENV_DIR="$JAWS_INSTALL_DIR/site"

echo "Installing into ${JAWS_CONFIG_DIR}"
echo "Creating paths and setting permissions"
FOLDERS="
JAWS_BIN_DIR
JAWS_CONFIG_DIR
JAWS_CROMWELL_EXECUTIONS_DIR
JAWS_DOWNLOADS_DIR
JAWS_INPUTS_DIR
JAWS_INSTALL_DIR
JAWS_LOGS_DIR
JAWS_SCRATCH_DIR
JAWS_SUPERVISOR_DIR
"
setup_dirs "$FOLDERS" "$JAWS_GROUP" 750

# these folders require special permissions:
chgrp "$JAWS_USER_GROUP" "$JAWS_INPUTS_DIR" && chmod 770 "$JAWS_INPUTS_DIR"
chgrp "$JAWS_USER_GROUP" "$JAWS_DOWNLOADS_DIR" && chmod 750 "$JAWS_DOWNLOADS_DIR"
chgrp "$JAWS_USER_GROUP" "$JAWS_CROMWELL_EXECUTIONS_DIR" && chmod 750 "$JAWS_CROMWELL_EXECUTIONS_DIR"

# check supervisord
[[ -n "$JAWS_LOAD_PYTHON" ]] && $JAWS_LOAD_PYTHON
if [[ ! -d "$JAWS_SUPERVISOR_DIR/bin" ]]; then
    echo "Installing supervisor"
    $JAWS_PYTHON -m venv "$JAWS_SUPERVISOR_DIR" && \
      . "$JAWS_SUPERVISOR_DIR/bin/activate" && \
      pip install --upgrade pip && \
      pip install supervisor && \
      deactivate
else
    echo "Stopping services"
    $JAWS_BIN_DIR/supervisorctl stop "jaws-site"
fi

echo "Generating virtual environment"
rm -rf ./*/dist/*
test -d "$JAWS_VENV_DIR" && rm -rf "$JAWS_VENV_DIR"
[[ -n "$JAWS_LOAD_PYTHON" ]] && $JAWS_LOAD_PYTHON
make pkg
$JAWS_PYTHON -m venv "$JAWS_VENV_DIR" && \
  . "$JAWS_VENV_DIR/bin/activate" && \
  pip install --upgrade pip && \
  pip install wheel && \
  pip install rpc/dist/* && \
  pip install site/dist/* && \
  deactivate

echo "Writing config files"
envsubst < "./test/integration/templates/site.conf" > "$JAWS_CONFIG_DIR/jaws-site.conf"
envsubst < "./test/integration/templates/supervisor.conf" > "$JAWS_CONFIG_DIR/supervisor.conf"
envsubst < "./test/integration/templates/supervisor.site.conf" > "$JAWS_CONFIG_DIR/supervisor.site.conf"
chmod 600 $JAWS_CONFIG_DIR/*

echo "Writing shims"
envsubst < "./test/integration/templates/rpc-server.sh" > "$JAWS_BIN_DIR/rpc-server"
envsubst < "./test/integration/templates/runs.sh" > "$JAWS_BIN_DIR/runs"
envsubst < "./test/integration/templates/transfers.sh" > "$JAWS_BIN_DIR/transfers"
envsubst < "./test/integration/templates/perf-metrics.sh" > "$JAWS_BIN_DIR/perf-metrics"
envsubst < "./test/integration/templates/supervisord.sh" > "$JAWS_BIN_DIR/supervisord"
envsubst < "./test/integration/templates/supervisorctl.sh" > "$JAWS_BIN_DIR/supervisorctl"
chmod 700 $JAWS_BIN_DIR/*

echo "Starting services"
$JAWS_BIN_DIR/supervisorctl start "jaws-site"

echo "END deploy-jaws"