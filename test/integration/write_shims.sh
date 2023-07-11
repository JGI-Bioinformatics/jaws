#!/usr/bin/env bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SERVICES=("rpc-server" "run-daemon" "transfer-daemon" "perf-metrics-daemon" "task-log")


function write_apptainer_run_shims {
  container_templ="$DIR/templates/container_runtime_templates/apptainer-run.sh"
  for service in "${SERVICES[@]}"; do
    export SERVICE="$service"
    envsubst < "$container_templ" > "$JAWS_BIN_DIR/$SERVICE"
  done
}

function write_apptainer_instance_shims {
  container_templ="$DIR/templates/container_runtime_templates/apptainer-instance.sh"
  for service in "${SERVICES[@]}"; do
    export SERVICE="$service"
    envsubst < "$container_templ" > "$JAWS_BIN_DIR/$SERVICE"
  done
}

function write_venv_shims {
    envsubst < "$DIR/templates/rpc-server.sh" > "$JAWS_BIN_DIR/rpc-server"
    envsubst < "$DIR/templates/task-log.sh" > "$JAWS_BIN_DIR/task-log"
    envsubst < "$DIR/templates/runs.sh" > "$JAWS_BIN_DIR/run-daemon"
    envsubst < "$DIR/templates/transfers.sh" > "$JAWS_BIN_DIR/transfer-daemon"
    envsubst < "$DIR/templates/perf-metrics.sh" > "$JAWS_BIN_DIR/perf-metrics-daemon"
    envsubst < "$DIR/templates/queue-wait.sh" > "$JAWS_BIN_DIR/queue-wait.sh"
}

function write_perlmutter_shims {
  envsubst < "$DIR/templates/jaws-site.sh" > "$JAWS_BIN_DIR/jaws-site-cronjob"
  envsubst < "$DIR/templates/jaws-perlmutter-gitlab-runner.sh" > "$JAWS_BIN_DIR/jaws-perlmutter-gitlab-runner"
}

function write_shims {
  local install_method="$1"
  case "$install_method" in
    "venv")
      write_venv_shims
      ;;
    "apptainer-run")
      write_apptainer_run_shims
      ;;
    "apptainer-instance")
      write_apptainer_instance_shims
      ;;
    "perlmutter")
      write_perlmutter_shims
      ;;
    *)
      echo "Unknown install method: $install_method"
      exit 1
      ;;
  esac
}