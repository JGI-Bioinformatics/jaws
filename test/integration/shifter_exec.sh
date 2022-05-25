#!/usr/bin/env bash
echo "Executing on $HOSTNAME" 1>&2
IMG=$1          # jfroula/bbtools@sha256:
SHELL=$2 # /bin/bash
SCRIPT=$3       # script
FAST_SCRATCH=$4
BIG_SCRATCH=$5

if [ $# -lt 5 ]; then
    echo "One or more arguments are missing."
    echo "Usage: $0 <docker image either tagged or sha@256> <path to bash> <script> <path to fast scratch> <path to big scratch>"
    exit 1
fi

MOUNT_FAST_SCRATCH=""
if [ ! -z $FAST_SCRATCH ] && [ -d $FAST_SCRATCH ]; then
    MOUNT_FAST_SCRATCH = "-V $FAST_SCRATCH:/fast_scratch"
fi

MOUNT_BIG_SCRATCH=""
if [ ! -z $BIG_SCRATCH ] && [ -d $BIG_SCRATCH ] ; then
    MOUNT_BIG_SCRATCH = "-V $BIG_SCRATCH:/big_scratch"
fi

IMG=${1}
REPO=$(echo $IMG | sed 's/@.*//')
HASH=$(echo $IMG | sed 's/.*@//')
ID=$(echo $IMG | sed 's/.*://')

# if using sha
#   IMG: jfroula/test@sha256:ef70f44e4d7cc28d40ff6117922583642919403dfac44d0a523f49c5f9b8993a
#   REPO: jfroula/test
#   HASH: sha256:ef70f44e4d7cc28d40ff6117922583642919403dfac44d0a523f49c5f9b8993a
#   ID: ef70f44e4d7cc28d40ff6117922583642919403dfac44d0a523f49c5f9b8993a
#
# if not using sha
#   IMG: jfroula/test:0.1.5
#   REPO: jfroula/test:0.1.5
#   HASH: jfroula/test:0.1.5
#   ID: 0.1.5

if [[ $HASH =~ "sha256" ]]; then
    shifter --image=id:${ID} echo testing to see if we already have image > /dev/null 2>&1
else
    shifter --image=${IMG} echo testing to see if we already have image > /dev/null 2>&1
fi
if [[ $? == 0 ]]; then
    echo "image already pulled: $IMG"
fi

# Run container script and catch exit code
if [[ $HASH =~ "sha256" ]]; then
    shifter --image=id:$ID $MOUNT_FAST_SCRATCH $MOUNT_BIG_SCRATCH $SHELL $SCRIPT
else
    shifter --image=$REPO $MOUNT_FAST_SCRATCH $MOUNT_BIG_SCRATCH $SHELL $SCRIPT
fi

export EXIT_CODE=$?

# Return with container exit code
exit $EXIT_CODE
