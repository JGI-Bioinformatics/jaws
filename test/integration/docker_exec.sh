#!/usr/bin/env bash
CWD=$1
DOCKER_CWD=$2
IMAGE=$3
JOB_SHELL=$4
SCRIPT=$5
FAST_SCRATCH=$6
BIG_SCRATCH=$7
echo "Executing on $HOSTNAME" 1>&2
echo "Begin execution at" `date +"%Y-%m-%d %H:%M:%S" -u` UTC 1>&2

MOUNT_FAST_SCRATCH=""
if [ ! -z $FAST_SCRATCH ] && [ -d $FAST_SCRATCH ]; then
    MOUNT_FAST_SCRATCH="--mount $FAST_SCRATCH:/fast_scratch"
fi

MOUNT_BIG_SCRATCH=""
if [ ! -z $BIG_SCRATCH ] && [ -d $BIG_SCRATCH ] ; then
    MOUNT_BIG_SCRATCH="--mount $BIG_SCRATCH:/big_scratch"
fi

# Run container script and catch exit code
docker run --rm -it $MOUNT_FAST_SCRATCH $MOUNT_BIG_SCRATCH -v $CWD:$DOCKER_CWD --workdir=$DOCKER_CWD --entrypoint=$JOB_SHELL $IMAGE $SCRIPT
EXIT_CODE=$?

echo "End execution at" `date +"%Y-%m-%d %H:%M:%S" -u` UTC 1>&2

# Return with container exit code
exit $EXIT_CODE
