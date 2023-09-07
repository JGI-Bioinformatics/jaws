#!/bin/bash -l

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export PYTHONIOENCODING=utf-8

source "$JAWS_VENV_DIR/bin/activate"
exec jaws-site --log "$JAWS_LOGS_DIR/task-logger-receive.log" --config "$JAWS_CONFIG_DIR/jaws-site.conf" --log-level $JAWS_LOG_LEVEL task-logger-receive
