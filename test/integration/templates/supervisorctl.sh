#!/bin/bash -l

source "$JAWS_SUPERVISOR_DIR/bin/activate"
exec supervisorctl -c "$JAWS_CONFIG_DIR/supervisor.conf" $@
