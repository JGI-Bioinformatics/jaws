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

function pullImage(){
    IMG=$1
    REPO=$2
    HASH=$3
    TAG=
    IMAGE=
    if [[ $HASH =~ "sha256" ]]; then
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
    
        IMAGE=${REPO}:${TAG}
    else
        IMAGE=${REPO}
    fi

    # Pull image by tag
    shifterimg pull ${IMAGE} > /dev/null 2>&1
    if [[ $? > 0 ]]; then
        echo "Invalid container name or failed to pull container, ${IMAGE}"
        exit 1
    else
        echo "successfully pulled image ${IMAGE}"
    fi
}

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
else
    pullImage $IMG $REPO $HASH
fi


# Check if monitoring script is installed
if [ -x "$(command -v pagurus)" ]; then
    # Start monitoring running, get monitor pid and sleep
    pagurus -u $USER -o $PWD/stats.csv &
    PID=$!
    sleep 2
fi

# Run container script and catch exit code
if [[ $HASH =~ "sha256" ]]; then
    shifter --image=id:$ID -V $2:$3 $4 $5
else
    shifter --image=$REPO -V $2:$3 $4 $5
fi

export EXIT_CODE=$?

# If PID is set then kill it
if [ -n "$PID" ]; then
    # Kill pagurus monitoring and wait for it to write stats.csv
    kill $PID
    sleep 2
fi

# Return with container exit code
exit $EXIT_CODE
