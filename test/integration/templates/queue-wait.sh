#!/usr/bin/env bash

#export LC_ALL=en_US.UTF-8
#export LANG=en_US.UTF-8
#export PYTHONIOENCODING=utf-8

small=${DOLLAR}($JAWS_QUEUE_WAIT_SMALL)
medium=${DOLLAR}($JAWS_QUEUE_WAIT_MEDIUM)
large=${DOLLAR}($JAWS_QUEUE_WAIT_LARGE)
xlarge=${DOLLAR}($JAWS_QUEUE_WAIT_XLARGE)
echo {\"small\":\"${DOLLAR}small\", \"medium\":\"${DOLLAR}medium\", \"large\":\"${DOLLAR}large\", \"xlarge\":\"${DOLLAR}xlarge\" }
