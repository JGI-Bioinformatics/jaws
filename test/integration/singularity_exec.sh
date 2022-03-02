#!/usr/bin/env bash
REMOTE_REF=$1
LOCAL_REF=$2
CWD=$3
DOCKER_CWD=$4
IMAGE=$5
JOB_SHELL=$6
SCRIPT=$7

if [ "$#" -ne 7 ]; then
    echo "Missing some arguments. There should be 7 arguments"
    echo "Usage: $0 <local refdata> <container refdata> <local cwd> <container cwd(same)> <image> <job shell> <script>"
    exit 1
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

# Check if monitoring script is installed
if [ -x "$(command -v pagurus)" ]; then
    # Start monitoring running, get monitor pid and sleep
    pagurus -u $USER -o $PWD/stats.csv &
    PID=$!
    sleep 2
fi

# Run container script and catch exit code
singularity exec --bind $REMOTE_REF:$LOCAL_REF --bind $CWD:$DOCKER_CWD $IMAGE $JOB_SHELL $SCRIPT
export EXIT_CODE=$?

# If PID is set then kill it
if [ -n "$PID" ]; then
    # Kill pagurus monitoring and wait for it to write stats.csv
    kill $PID
    sleep 2
fi

# Return with container exit code
exit $EXIT_CODE
