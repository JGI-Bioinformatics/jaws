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
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export SITE_NAME=$SITE_NAME
export DEPLOY_NAME=$DEPLOY_NAME
export DEPLOY_DIR=$DEPLOY_DIR

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
    cmd="$SRC_DIR/setup_jaws $SITE_NAME $DEPLOY_NAME $INSTALL_DIR"
else
    cmd="$SRC_DIR/setup_jaws $SITE_NAME $DEPLOY_NAME $INSTALL_DIR -c $4"
fi
echo
echo "Setting up deployment environment ..."
echo "$cmd"
eval $cmd

echo
echo "Running jaws services..."
echo "$DEPLOY_DIR/run_jaws"
$DEPLOY_DIR/run_jaws
