#!/bin/bash -l

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export PYTHONIOENCODING=utf-8

source "$JAWS_VENV_DIR/bin/activate"
exec jaws-site --log "$JAWS_LOGS_DIR/site-perf-metrics-daemon.log" --config "$JAWS_CONFIG_DIR/jaws-site.conf" --log-level $JAWS_LOG_LEVEL perf-metrics-daemon
