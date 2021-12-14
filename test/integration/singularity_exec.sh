#!/usr/bin/env bash
REMOTE_REF=$1
LOCAL_REF=$2
CWD=$3
DOCKER_CWD=$4
IMAGE=$5
JOB_SHELL=$6
SCRIPT=$7
echo "Executing on $HOSTNAME" 1>&2
singularity exec --bind $REMOTE_REF:$LOCAL_REF --bind $CWD:$DOCKER_CWD $IMAGE $JOB_SHELL $SCRIPT
exit $?
