#!/bin/csh -f
echo -n "hostname="; hostname;
echo "pid=$$"
sleep 10
ps -ef | grep tmfq-worker
