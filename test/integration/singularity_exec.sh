#!/usr/bin/env bash
CWD=$1
DOCKER_CWD=$2
IMAGE=$3
JOB_SHELL=$4
SCRIPT=$5
FAST_SCRATCH=$6
BIG_SCRATCH=$7

if [ "$#" -ne 7 ]; then
    echo "Missing some arguments. There should be 7 arguments"
    echo "Usage: $0 <local cwd> <container cwd(same)> <image> <job shell> <script> <fast_scratch> <big_scratch>"
    exit 1
fi

MOUNT_FAST_SCRATCH=""
if [ ! -z $FAST_SCRATCH ] && [ -d $FAST_SCRATCH ]; then
    MOUNT_FAST_SCRATCH = "--bind $FAST_SCRATCH:/fast_scratch"
fi

MOUNT_BIG_SCRATCH=""
if [ ! -z $BIG_SCRATCH ] && [ -d $BIG_SCRATCH ] ; then
    MOUNT_BIG_SCRATCH = "--bind $BIG_SCRATCH:/big_scratch"
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
singularity exec $MOUNT_FAST_SCRATCH $MOUNT_BIG_SCRATCH --bind $CWD:$DOCKER_CWD $IMAGE $JOB_SHELL $SCRIPT
export EXIT_CODE=$?

# Return with container exit code
exit $EXIT_CODE
