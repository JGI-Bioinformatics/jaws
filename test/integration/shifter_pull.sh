#!/usr/bin/env bash
# This script will:
# 1) check to see if it is cached already.
# 2) if not, figure out the right version to pull

IMG=$1          # jfroula/bbtools@sha256:

if [ ! $IMG ]; then
    echo "Usage: $0 <docker image either versioned tagged or with sha@256>"
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

set -eo pipefail
echo IMG: $IMG
echo REPO: $REPO
echo HASH: $HASH
echo ID: $ID
img_exists=$(shifter --image=id:${ID} echo imageExists)
echo img_exists: $img_exists
if [[ $img_exists -eq "imageExists" ]]; then
    echo "image already pulled: $IMG"
    exit 0
fi

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
shifterimg pull ${IMAGE}
