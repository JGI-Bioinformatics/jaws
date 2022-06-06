#!/bin/bash

echo "Stopping worker $(hostname)!"
echo $SLURM_HOST

# Kill condor_master and let it kill all other condor processes
kill $(ps aux | grep -v grep | grep -i condor_master | awk '{print $2}')
# Kills pagurus
kill $(ps aux | grep -v grep | grep -i pagurus | awk '{print $2}')

sleep 10