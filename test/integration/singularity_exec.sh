#!/usr/bin/env bash
REMOTE_REF=$1
LOCAL_REF=$2
CWD=$3
DOCKER_CWD=$4
IMAGE=$5
JOB_SHELL=$6
SCRIPT=$7

## VERIFY REQUIRED VARS ARE DEFINED
# Exits if any required var is undefined.
REQUIRED_VARS="
SINGULARITY_CACHEDIR
SINGULARITY_PULLFOLDER
SINGULARITY_TMPDIR
SINGULARITY_LOCALCACHEDIR
FLOCK_DIR
"
RESULT=0
for VAR in $REQUIRED_VARS; do
  if [[ -z ${!VAR} ]]; then
    RESULT=1
    echo "Missing env var: $VAR">&2
  fi
done
[ $RESULT -eq 0 ] || exit 1

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
