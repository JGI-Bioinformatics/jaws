#!/usr/bin/env bash
LOCAL_REF=$1
DOCKER_REF=$2
CWD=$3
DOCKER_CWD=$4
IMAGE=$5
JOB_SHELL=$6
SCRIPT=$7

if [ "$#" -ne 7 ]; then
    echo "Missing some arguments. There should be 7 arguments"
    echo "Usage: $0 <local refdata> <container refdata> <local cwd> <container cwd(same)> <image> <job shell> <script>"
    exit 1
fi

if [[ ! -d $CWD ]]; then
    echo "The CWD dir: $CWD is not a valid directory"
    exit 1
fi

if [[ ! -d $LOCAL_REF ]]; then
    echo "The Local reference dir: $LOCAL_REF is not a valid directory"
    exit 1
fi

if [[ ! -e $SCRIPT ]]; then
    echo "The SCRIPT: $SCRIPT doesn't exist"
    exit 1
fi

echo "Executing on $HOSTNAME" 1>&2

# Run container script and catch exit code
singularity exec --bind $LOCAL_REF:$DOCKER_REF --bind $CWD:$DOCKER_CWD $IMAGE $JOB_SHELL $SCRIPT
export EXIT_CODE=$?


# Return with container exit code
exit $EXIT_CODE
