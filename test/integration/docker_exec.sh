#!/usr/bin/env bash
REMOTE_REF=$1
LOCAL_REF=$2
CWD=$3
DOCKER_CWD=$4
IMAGE=$5
JOB_SHELL=$6
SCRIPT=$7
echo "Executing on $HOSTNAME" 1>&2

# Check if monitoring script is installed
if [ -x "$(command -v pagurus)" ]; then
    # Start monitoring running, get monitor pid and sleep
    pagurus -u $USER -o $PWD/stats.csv &
    PID=$!
    sleep 2
fi

# Run container script and catch exit code
docker run --rm -it -v $REMOTE_REF:$LOCAL_REF -v $CWD:$DOCKER_CWD --workdir=$DOCKER_CWD --entrypoint=$JOB_SHELL $IMAGE $SCRIPT
EXIT_CODE=$?

# If PID is set then kill it
if [ -n "$PID" ]; then
    # Kill pagurus monitoring and wait for it to write stats.csv
    kill $PID
    sleep 2
fi

# Return with container exit code
exit $EXIT_CODE
