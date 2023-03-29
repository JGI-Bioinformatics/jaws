#!/usr/bin/env bash -l

#export LC_ALL=en_US.UTF-8
#export LANG=en_US.UTF-8
#export PYTHONIOENCODING=utf-8
#
#source "$JAWS_VENV_DIR/bin/activate"
#exec jaws-site --log "$JAWS_LOGS_DIR/s

medium=$($JAWS_QUEUE_WAIT_MEDIUM)
xlarge=$($JAWS_QUEUE_WAIT_XLARGE)
echo {\"medium:\" $medium, \"xlarge:\" $xlarge}
