#!/usr/bin/env bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) # noqa

function write_jaws_configs {
  echo "Writing config files"
  envsubst < "$DIR/templates/site.conf" > "$JAWS_CONFIG_DIR/jaws-site.conf"
  envsubst < "$DIR/templates/site.env.templ" > "$JAWS_CONFIG_DIR/site.env"
  envsubst < "$DIR/templates/supervisor.conf" > "$JAWS_CONFIG_DIR/supervisor.conf"
  chmod 600 $JAWS_CONFIG_DIR/*.conf

  envsubst < "$DIR/templates/supervisord.sh" > "$JAWS_BIN_DIR/supervisord"
  envsubst < "$DIR/templates/supervisorctl.sh" > "$JAWS_BIN_DIR/supervisorctl"
  chmod 700 $JAWS_BIN_DIR/*

  if [[ "$JAWS_SETFACL" -eq 1 ]]; then
      echo "Setup file ACL rules"
      FACL_SCRIPT="$JAWS_SCRATCH_DIR/setup_facl.sh"
      envsubst < "$DIR/templates/setup_facl.sh" > "$FACL_SCRIPT"
      chmod 700 "$FACL_SCRIPT"
      "$FACL_SCRIPT"
  fi
}

function write_supervisor_configs {
  envsubst < "$DIR/templates/supervisor.site.conf" > "$JAWS_CONFIG_DIR/supervisor.site.conf"
  chmod 600 $JAWS_CONFIG_DIR/*.conf
  chmod 700 $JAWS_BIN_DIR/*
}

function write_systemd_configs {
  echo "Writting systemd service files"
  SERVICES=("rpc-server" "run-daemon" "transfer-daemon" "perf-metrics-daemon" "task-logger-receive")
  service_dir="${HOME}/.config/systemd/user"
  for service in "${SERVICES[@]}"; do
    export SERVICE="$service"
    envsubst < "$DIR/templates/jaws-site.service.templ" > "${service_dir}/${SERVICE}-${JAWS_DEPLOYMENT_NAME}.service"
  done
}

function write_configs {
  local launch_tool="$1"

  # always write the jaws configs
  write_jaws_configs

  case "$launch_tool" in
    "supervisor")
      write_supervisor_configs
      ;;
    "systemd")
      write_systemd_configs
      ;;
    "scrontab")
      write_supervisor_configs
      ;;
    *)
      echo "Unknown launch tool: $launch_tool"
      exit 1
      ;;
  esac

}
