#!/usr/bin/env bash -l

#export LC_ALL=en_US.UTF-8
#export LANG=en_US.UTF-8
#export PYTHONIOENCODING=utf-8

small=$($JAWS_QUEUE_WAIT_SMALL)
medium=$($JAWS_QUEUE_WAIT_MEDIUM)
large=$($JAWS_QUEUE_WAIT_LARGE)
xlarge=$($JAWS_QUEUE_WAIT_XLARGE)
echo {\"small\":\"${DOLLAR}small\", \"medium\":\"${DOLLAR}medium\", \"large\":\"${DOLLAR}large\", \"xlarge\":\"${DOLLAR}xlarge\" }
