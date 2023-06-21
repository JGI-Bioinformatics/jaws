#!/usr/bin/env bash

set -euxo pipefail

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$DIR/utils.sh"
source "$DIR/write_configs.sh"
source "$DIR/write_shims.sh"
source "$DIR/install_jaws_site.sh"

# Simple function that accepts the arguments for the launch tool (eg supervisor, systemd, etc) and the
# install method (eg venv, apptainer-sif, apptainer-docker) and calls the functions  loaded from sourced files
# to write the configs, write the shims, and install the application.
function deploy {
  launch_tool="$1"
  install_method="$2"
  echo "LAUNCH_TOOL=$launch_tool"
  write_configs
  write_shims "$launch_tool" "$install_method"
  install_jaws_site "$install_method"
  "${DIR}/deploy_tools/${launch_tool}.sh" start
}

deploy "$@"