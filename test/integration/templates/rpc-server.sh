#!/bin/bash -l

set -a
source "$JAWS_CONFIG_DIR/site.env"
set +a

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export PYTHONIOENCODING=utf-8

source "$JAWS_VENV_DIR/bin/activate"
export "$JAWS_BIN_DIR":$DOLLAR{PATH}
exec jaws-site --log "$JAWS_LOGS_DIR/site-rpc.log" --config "$JAWS_CONFIG_DIR/jaws-site.conf" --log-level $JAWS_LOG_LEVEL rpc-server
