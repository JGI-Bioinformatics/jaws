#!/bin/bash

# This script is only used when jaws-site is submitted to a cluster scheduler (e.g. SLURM)
# rather than being installed on a dedicated server.

JAWS_BIN_DIR="$JAWS_BIN_DIR"
JAWS_LOGS_DIR="$JAWS_LOGS_DIR"
JAWS_GITLAB_RUNNER="$JAWS_GITLAB_RUNNER"
JAWS_GITLAB_RUNNER_CONFIG="$JAWS_GITLAB_RUNNER_CONFIG"
LOGFILE="$JAWS_LOGS_DIR/jaws-site.log"


echo "Starting jaws-site on `hostname -i`" > "${DOLLAR}LOGFILE"


# shutdown cleanly
function stopall {
    echo "Stopping services" > "${DOLLAR}LOGFILE"
    "${DOLLAR}JAWS_BIN_DIR/supervisorctl" stop all 2>&1 >> "${DOLLAR}LOGFILE"
    echo "Shutdown supervisor" >> "${DOLLAR}LOGFILE"
    "${DOLLAR}JAWS_BIN_DIR/supervisorctl" shutdown 2>&1 >> "${DOLLAR}LOGFILE"
    echo "Stopping gitlab-runner" >> "${DOLLAR}LOGFILE"
    "${DOLLAR}JAWS_GITLAB_RUNNER" stop -c "${DOLLAR}JAWS_GITLAB_RUNNER_CONFIG" 2>&1 >> "${DOLLAR}LOGFILE"
    echo "Done" >> "${DOLLAR}LOGFILE"
}
trap stopall EXIT


# Start gitlab-runner
# This is not part of the CI/CD and must be setup manually first.
"${DOLLAR}JAWS_GITLAB_RUNNER" run -c "${DOLLAR}JAWS_GITLAB_RUNNER_CONFIG"

# Start supervisord
"${DOLLAR}JAWS_BIN_DIR/supervisord"