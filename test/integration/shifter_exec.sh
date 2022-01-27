#!/usr/bin/env bash
echo "Executing on $HOSTNAME" 1>&2
IMG=$1          # jfroula/bbtools@sha256:
DB=$2           # /global/dna/shared/databases/jaws/refdata/
MOUNT=$3        # /refdata
SHELL_SCRIPT=$4 # /bin/bash
SCRIPT=$5       # script

if [ $# -lt 5 ]; then
    echo "One or more arguments are missing."
    echo "Usage: $0 <docker image either tagged or sha@256> <path to refdata on nfs> <path to refdata in container> <path to bash> <script>"
    exit 1
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
#   Tue Jan 18 14:05:51 jaws@cori20 /global/cfs/cdirs/jaws/jaws-install/jaws-prod$ ./shifter_pull.sh jfroula/test:0.1.5
#   IMG: jfroula/test:0.1.5
#   REPO: jfroula/test:0.1.5
#   HASH: jfroula/test:0.1.5
#   ID: 0.1.5

if [[ $HASH =~ "sha256" ]]; then
    ID=id:$ID
else
    ID=$1
fi

# Check if monitoring script is installed
if [ -x "$(command -v pagurus)" ]; then
    # Start monitoring running, get monitor pid and sleep
    pagurus -u $USER -o $PWD/stats.csv &
    PID=$!
    sleep 2
fi

# Run container script and catch exit code
shifter --image=$ID -V $2:$3 $4 $5
export EXIT_CODE=$?

# If PID is set then kill it
if [ -n "$PID" ]; then
    # Kill pagurus monitoring and wait for it to write stats.csv
    kill $PID
    sleep 2
fi

# Return with container exit code
exit $EXIT_CODE
