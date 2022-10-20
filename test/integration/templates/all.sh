export JAWS_INSTALL_DIR="$JAWS_INSTALL_BASEDIR/${JAWS_SITE_NAME}-${JAWS_DEPLOYMENT_NAME}"
export JAWS_CONFIG_DIR="$JAWS_INSTALL_DIR/config"
export JAWS_BIN_DIR="$JAWS_INSTALL_DIR/bin"
export JAWS_LOGS_DIR="$JAWS_INSTALL_DIR/log"
export JAWS_SUPERVISOR_DIR="$JAWS_INSTALL_DIR/supervisor"
if [[ $JAWS_SCRATCH_BASEDIR =~ ^s3:\/\/ ]]; then
    # EACH AWS DEPLOYMENT HAS A SEPARATE S3 BUCKET
    export JAWS_SCRATCH_DIR="${JAWS_SCRATCH_BASEDIR}-${JAWS_DEPLOYMENT_NAME}"
else
    # EACH SERVER DEPLOYMENT HAS A SEPARATE NFS SUBDIR
    export JAWS_SCRATCH_DIR="$JAWS_SCRATCH_BASEDIR/jaws-${JAWS_DEPLOYMENT_NAME}"
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
export JAWS_GLOBUS_ROOT_DIR

# virtual envs
export JAWS_VENV_DIR="$JAWS_INSTALL_DIR/site"
