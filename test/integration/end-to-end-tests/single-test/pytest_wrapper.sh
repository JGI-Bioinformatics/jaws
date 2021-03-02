#!/bin/bash

TEST_FOLDER=$1

# check we have all arguments and environmental variables set
if [[ ! $TEST_FOLDER ]]; then
	echo "Usage: $0 <folder-with-tests>"
	exit 1
fi


set -euox pipefail

## VALIDATE INPUT VARS
## print list of all missing vars (instead of exiting on first undefined var)
## exit with 1 if any vars are missing; 0 otherwise
REQUIRED_VARS="
JAWS_SITE
DEPLOYMENT_NAME
JAWS_TEST_TOKEN
"
RESULT=0
echo "VALIDATING REQUIRED VARS"
for VAR in $REQUIRED_VARS; do
  if [[ -z ${!VAR} ]]; then
    echo "Missing env var: $VAR">&2
    RESULT=1
  fi
  echo "... $VAR : OK"
done
[ $RESULT -eq 0 ] || exit $RESULT


## SITE-SPECIFIC VARIABLES:
## This will define global vars with a SITE_ prefix, with values from variables referenced by the JAWS_SITE
## variable.  For instance, if JAWS_SITE="JGI", setting GLOBUS_EP site var will define a SITE_GLOBUS_EP var
## with the value from JGI_GLOBUS_EP.
## This allows the site installation scripts to use the SITE_* vars.
echo "DEFINING SITE VARS"
function set_site_var {
  VAR_NAME="$1"
  SRC_VAR_NAME="${JAWS_SITE}_${VAR_NAME}"
  if [ -z ${!SRC_VAR_NAME+xxx} ]; then
    echo "Missing env var: $SRC_VAR_NAME">&2
    exit 1
  fi
  SITE_VAR_NAME="SITE_${VAR_NAME}"
  VALUE="${!SRC_VAR_NAME}"
  echo "... $SITE_VAR_NAME"
  export $SITE_VAR_NAME
  printf -v "$SITE_VAR_NAME" "%s" "$VALUE"
}
SITE_VARS="
INSTALL_BASEDIR
JAWS_SW_BASEDIR
"
for SITE_VAR in $SITE_VARS; do
  set_site_var $SITE_VAR
done


## DEFINE VARS FROM OTHER VARS
echo "DEFINING PATHS"
export INSTALL_DIR="$SITE_INSTALL_BASEDIR/jaws-$DEPLOYMENT_NAME"
export SITE_CLIENT_INSTALL_DIR="$SITE_JAWS_SW_BASEDIR/jaws-$DEPLOYMENT_NAME"

# source venv
source "$SITE_CLIENT_INSTALL_DIR/bin/activate"

# write token to homedir
echo << END > ~/jaws.conf
{
	token: $JAWS_TEST_TOKEN
}
END 

chmod 600 ~/jaws.conf

cd test/integration/end-to-end-tests
pytest --verbose --env $DEPLOYMENT_NAME --site $JAWS_SITE $TEST_FOLDER
