#!/usr/bin/env bash
CWD=$1
DOCKER_CWD=$2
IMAGE=$3
JOB_SHELL=$4
SCRIPT=$5
REF_DATA_DIR=$6
FAST_SCRATCH=$7
BIG_SCRATCH=$8

if [ "$#" -ne 8 ]; then
    echo "Missing some arguments. There should be 8 arguments"
    echo "Usage: $0 <local cwd> <container cwd(same)> <image> <job shell> <script> <ref_data> <fast_scratch> <big_scratch>"
    exit 1
fi

MOUNT_REF_DATA=""
if [ ! -z $REF_DATA_DIR ] && [ -d $REF_DATA_DIR ]; then
    MOUNT_REF_DATA="--bind $REF_DATA_DIR:/refdata"
fi

MOUNT_FAST_SCRATCH=""
if [ ! -z $FAST_SCRATCH ] && [ -d $FAST_SCRATCH ]; then
    MOUNT_FAST_SCRATCH="--bind $FAST_SCRATCH:/fast_scratch"
fi

MOUNT_BIG_SCRATCH=""
if [ ! -z $BIG_SCRATCH ] && [ -d $BIG_SCRATCH ] ; then
    MOUNT_BIG_SCRATCH="--bind $BIG_SCRATCH:/big_scratch"
fi

if [[ ! -d $CWD ]]; then
    echo "The CWD dir: $CWD is not a valid directory"
    exit 1
fi

if [[ ! -e $SCRIPT ]]; then
    echo "The SCRIPT: $SCRIPT doesn't exist"
    exit 1
fi

echo "Executing on $HOSTNAME" 1>&2
echo "Begin execution at" `date +"%Y-%m-%d %H:%M:%S" -u` UTC 1>&2

# Run container script and catch exit code
singularity exec $MOUNT_REF_DATA $MOUNT_FAST_SCRATCH $MOUNT_BIG_SCRATCH --bind $CWD:$DOCKER_CWD $IMAGE $JOB_SHELL $SCRIPT
export EXIT_CODE=$?

echo "End execution at" `date +"%Y-%m-%d %H:%M:%S" -u` UTC 1>&2

# Return with container exit code
exit $EXIT_CODE
