#!/bin/bash

set -eo pipefail

# INPUTS:
# $1 = site name (e.g. cori, jgi)
# $2 = deploy_env name (e.g. dev, staging, prod)
# $3 = install_dir
# $4 = config file (deploy json file - optional)

SITE_NAME=$1
DEPLOY_NAME=$2
INSTALL_DIR=$(readlink -f $3)
DEPLOY_DIR="$INSTALL_DIR/deploy"
SRC_DIR=$(pwd)

export SITE_NAME=$SITE_NAME
export DEPLOY_NAME=$DEPLOY_NAME
export DEPLOY_DIR=$DEPLOY_DIR
export SRC_DIR=$SRC_DIR

# ** FOR TESTING ONLY **
export JAWS_VERSION="2.0.1"
export JAWS_CENTRAL_DB_PW="MBQG0!7jAxk1M2R6Ufh"
export JAWS_CENTRAL_RMQ_PW="vG3a8KY5T36Vs6FrYJ8"
export JAWS_SITE_CORI_DB_PW="!DG3Cw7RtaYkdVrd8hv"
export JAWS_SITE_CORI_RMQ_PW="vG3a8KY5T36Vs6FrYJ8"
export JAWS_CROMWELL_CORI_DB_PW="!DG3Cw7RtaYkdVrd8hv"
export JAWS_JTM_CORI_DB_PW="!DG3Cw7RtaYkdVrd8hv"
export JAWS_JTM_RMQ_CORI_PW="vG3a8KY5T36Vs6FrYJ8"
export JAWS_SUPERVISORD_PW="zJcKZSxgPYSNF3PD"
# END TESTING

test -d "$INSTALL_DIR" && rm -rf "$INSTALL_DIR"
mkdir "$INSTALL_DIR"
mkdir "$DEPLOY_DIR"

echo "Creating python env for deployment script ..."
module load python  # NEED TO MOVE THIS TO .gitlab-ci.yml
python3 -m venv $DEPLOY_DIR/venv && \
  . "$DEPLOY_DIR/venv/bin/activate" && \
  pip install --upgrade pip && \
  pip install jinja2 wheel

if [[ -z "$4" ]];then
    cmd="./test/integration/deploy/setup_jaws $SITE_NAME $DEPLOY_NAME $INSTALL_DIR"
else
    cmd="./test/integration/deploy/setup_jaws $SITE_NAME $DEPLOY_NAME $INSTALL_DIR -c $4"
fi
echo
echo "Setting up deployment environment ..."
echo "$cmd"
eval $cmd

echo
echo "Running jaws services..."
echo "$DEPLOY_DIR/run_jaws"
$DEPLOY_DIR/run_jaws