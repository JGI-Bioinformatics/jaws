#!/bin/bash

# This script is only used when jaws-site is submitted to a cluster scheduler (e.g. SLURM)
# rather than being installed on a dedicated server.

JAWS_BIN_DIR="$JAWS_BIN_DIR"
JAWS_LOGS_DIR="$JAWS_LOGS_DIR"
JAWS_GITLAB_RUNNER="$JAWS_GITLAB_RUNNER"
JAWS_GITLAB_RUNNER_CONFIG="$JAWS_GITLAB_RUNNER_CONFIG"
LOGFILE="${DOLLAR}JAWS_LOGS_DIR/jaws-site.log"


echo "Starting jaws-site on `hostname -i`" > "$LOGFILE"


# shutdown cleanly
function stopall {
    echo "Stopping services" > "$LOGFILE"
    "${DOLLAR}JAWS_BIN_DIR/supervisorctl" stop all 2>&1 >> "$LOGFILE"
    echo "Shutdown supervisor" >> "$LOGFILE"
    "${DOLLAR}JAWS_BIN_DIR/supervisorctl" shutdown 2>&1 >> "$LOGFILE"
    echo "Stopping gitlab-runner" >> "$LOGFILE"
    "${DOLLAR}JAWS_GITLAB_RUNNER" stop -c "${DOLLAR}JAWS_GITLAB_RUNNER_CONFIG" 2>&1 >> "$LOGFILE"
    echo "Done" >> "$LOGFILE"
}
trap stopall EXIT


# Start gitlab-runner
# This is not part of the CI/CD and must be setup manually first.
"${DOLLAR}JAWS_GITLAB_RUNNER" run -c "${DOLLAR}JAWS_GITLAB_RUNNER_CONFIG"

# Start supervisord
exec "${DOLLAR}JAWS_BIN_DIR/supervisord"
