#!/bin/bash

# This is to request a worklfow node from Perlmutter and start a gitlab-runner
# This gitlab-runner can be shared among all possible JAWS sites like "perlmutter" or "nmdc"
# This also can be shared among all deployments
JAWS_GITLAB_RUNNER="$JAWS_GITLAB_RUNNER"
JAWS_GITLAB_RUNNER_CONFIG="$JAWS_GITLAB_RUNNER_CONFIG"

# Start gitlab-runner
"${DOLLAR}JAWS_GITLAB_RUNNER" run -c "${DOLLAR}JAWS_GITLAB_RUNNER_CONFIG"