#!/usr/bin/env bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

function write_configs {
  echo "Writing config files"
  envsubst < "$DIR/templates/site.conf" > "$JAWS_CONFIG_DIR/jaws-site.conf"
  envsubst < "$DIR/templates/site.env.templ" > "$JAWS_CONFIG_DIR/site.env"
  envsubst < "$DIR/templates/supervisor.conf" > "$JAWS_CONFIG_DIR/supervisor.conf"
  envsubst < "$DIR/templates/supervisor.site.conf" > "$JAWS_CONFIG_DIR/supervisor.site.conf"
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