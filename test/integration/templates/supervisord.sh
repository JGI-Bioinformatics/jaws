#!/bin/bash -l

$JAWS_LOAD_PYTHON
source "$JAWS_SUPERVISOR_DIR/bin/activate"
exec supervisord -c "$JAWS_CONFIG_DIR/supervisor.conf" $@
