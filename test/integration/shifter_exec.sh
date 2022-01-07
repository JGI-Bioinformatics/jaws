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

function pull_by_id() {
    # This will:
    # 1) check to see if it is cached already.
    # 2) if not, figure out the right version to pull

    set -e
    IMG=${1}
    REPO=$(echo $IMG | sed 's/@.*//')
    HASH=$(echo $IMG | sed 's/.*@//')
    ID=$(echo $IMG | sed 's/.*://')

    # Try running it...
    set +e
    shifter --image=id:${ID} echo yes >/dev/null 2>&1 && return 0
    set -e

    #Try to figure out the version to pull
    RT=$(skopeo inspect docker://${IMG} | jq .RepoTags)
    for ttag in $(echo $RT | sed 's/[",[]//g'); do
        digest=$(skopeo inspect docker://${REPO}:${ttag} | jq .Digest | sed 's/"//g')
        if [ "$digest" == "$HASH" ]; then
            TAG=$ttag
            break
        fi
    done
    if [ -z $TAG ]; then
        echo "Unable to determine image version" 1>&2
        exit 1
    fi

    # Pull image by tag
    shifterimg pull ${REPO}:${TAG} 1>&2

    # Get the ID
    # All variables are global by default in bash so we can access the
    # ID variable outside the function now.
    ID=$(shifterimg lookup ${REPO}:${TAG})
}

if [ $(echo ${IMG} | grep -c sha256) -gt 0 ]; then
    pull_by_id ${IMG}
    ID=id:$ID
else
    ID=$1
fi

# Start monitoring running, get monitor pid and sleep
pagurus -u $USER -o $PWD/stats.csv &
export PID=$!
sleep 2

# Run container script and catch exit code
shifter --image=$ID -V $2:$3 $4 $5
export EXIT_CODE=$?

# Kill pagurus monitoring and wait for it to write stats.csv
kill $PID
sleep 5

# Return with container exit code
exit $EXIT_CODE
