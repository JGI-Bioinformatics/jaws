#!/bin/bash

# This script is only used when jaws-site is submitted to a cluster scheduler (e.g. SLURM)
# rather than being installed on a dedicated server.

# Start gitlab-runner
# This is not part of the CI/CD and must be setup manually first.
"$JAWS_GITLAB_RUNNER" run -c "$JAWS_GITLAB_RUNNER_CONFIG"

# Start supervisord
export JAWS_BIN_DIR="$JAWS_BIN_DIR"
exec ${DOLLAR}JAWS_BIN_DIR/supervisord
