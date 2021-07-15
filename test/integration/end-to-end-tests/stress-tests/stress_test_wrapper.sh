#!/bin/bash
# This is a wrapper to set the required variables so that we can source the appropriate 
# JAWS environmental file (i.e jaws-dev.sh) so that the pytest(s) that use jaws client will work.
# After the vars are all set, the pytest(s) are run.
#
# Assumptions:
# ------------
# It is assumed that the following variables are exported in the 
# shell environment. This is usually done through the CI/CD infrastructure.
# DEPLOYMENT_NAME
# JAWS_SITE
# CORI_BASEDIR
# JAWS_TEST_TOKEN
# 
# note that JAWS_SITE may be a list of sites (i.e. JAWS_SITE="cori jgi")

TEST_SCRIPT=$1

# check we have all arguments and environmental variables set
if [[ ! $TEST_SCRIPT ]]; then
	echo "Usage: $0 <pytest script>"
	exit 1
fi

## VALIDATE INPUT VARS
## print list of all missing vars (instead of exiting on first undefined var)
## exit with 1 if any vars are missing; 0 otherwise
REQUIRED_VARS="
DEPLOYMENT_NAME
JAWS_SITE
CORI_BASEDIR
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


# write token to jaws.conf
# [USER]
# token =
echo -e "[USER]\ntoken = $JAWS_TEST_TOKEN" > ~/jaws.conf
chmod 600 ~/jaws.conf

echo "source $CORI_BASEDIR/jaws-${DEPLOYMENT_NAME}/jaws-${DEPLOYMENT_NAME}.sh"
cd test/integration/end-to-end-tests
echo /global/cfs/projectdirs/jaws/jaws-dev/bin/pytest --verbose --env ${DEPLOYMENT_NAME} --site-list \"${JAWS_SITE}\" ${TEST_SCRIPT}
	 /global/cfs/projectdirs/jaws/jaws-dev/bin/pytest --verbose --env ${DEPLOYMENT_NAME} --site-list \"${JAWS_SITE}\" ${TEST_SCRIPT}

