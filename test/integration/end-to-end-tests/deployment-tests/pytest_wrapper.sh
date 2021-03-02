#!/bin/bash

ENV=$1
SITE=$2
TEST_FOLDER=$3

# check we have all arguments and environmental variables set
if [[ ! $ENV ]] || [[ ! $SITE ]] || [[ ! $TEST_FOLDER ]]; then
	echo "Usage: $0 <[prod|staging|dev]> <[cori|jgi]> <folder-with-tests>"
	exit 1
fi

if [[ ! -d $SITE_JAWS_SW_BASEDIR ]]; then
	echo "SITE_JAWS_SW_BASEDIR environmental variable not set"
	exit 1
fi
if [[ ! $DEPLOYMENT_NAME ]]; then
	echo "DEPLOYMENT_NAME environmental variable not set"
	exit 1
fi
if [[ ! $JAWS_TEST_USER_TOKEN ]]; then
	echo "JAWS_TEST_USER_TOKEN environmental variable not set"
	exit 1
fi


set -euox pipefail

# source venv
source $SITE_JAWS_SW_BASEDIR/jaws-$DEPLOYMENT_NAME/bin/activate

# write token to homedir
echo << END > ~/jaws.conf
{
	token: $JAWS_TEST_USER_TOKEN
}
END 

chmod 600 ~/jaws.conf

cd test/integration/end-to-end-tests
pytest --verbose --env $ENV --site $SITE $TEST_FOLDER
