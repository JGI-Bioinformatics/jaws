#!/bin/bash

set -e
#source ~/venv/bin/activate
# set env activation

#MODE="prod"
#MODE="dev"
MODE=$1
_now=$(date +"%Y_%m_%d_%I_%M_%p")
if [ -z "$JTM_HOST_NAME" ]; then echo "Please set JTM_HOST_NAME env variable"; exit -1; fi

#if [ -z "$NERSC_HOST" ]; then export NERSC_HOST=`cat /etc/clustername`; fi
#if [ -z "$HOSTNAME" ]; then export NERSC_HOST=`hostname`; fi
if [ -z "$SCRATCH" ]; then export $SCRATCH=/tmp; fi

if [  "$JTM_HOST_NAME" = "cori" ]; then
    # cori or cori20
    # if [  "$MODE" = "prod" ]; then
    #     ENV_ACTIVATION="source ~/venv/bin/activate"
    # else
    #     ENV_ACTIVATION="source ~/venv-dev/bin/activate"
    # fi
    # ENV_ACTIVATION="source /global/project/projectdirs/jaws_jtm/anaconda3/bin/activate /global/project/projectdirs/jaws_jtm/$MODE/jtm"
    ENV_ACTIVATION="source /global/cfs/projectdirs/jaws_jtm/anaconda3/bin/activate /global/cfs/projectdirs/jaws_jtm/$MODE/jtm"
elif [  "$JTM_HOST_NAME" = "jgi" ]; then
    # lrc-services.lbl.gov
    ENV_ACTIVATION="source /global/home/groups-sw/lr_jgicloud/anaconda3/bin/activate /global/home/groups-sw/lr_jgicloud/$MODE/jtm"
else
    # local
    if [  "$MODE" = "prod" ]; then
        ENV_ACTIVATION="source ~/venv/bin/activate"
    else
        ENV_ACTIVATION="source ~/venv-dev/bin/activate"
    fi
fi
eval "$ENV_ACTIVATION" && \
python <<CODE
import sys, os, glob
print sum([i+i**i for i in range(20000)])
CODE
#sleep 30

ls
