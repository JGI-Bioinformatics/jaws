#!/bin/bash

echo "Stopping worker $(hostname)!"
echo $SLURM_HOST

# Kill condor_master and let it kill all other condor processes
kill $(ps aux | grep -v grep | grep -i condor_master | awk '{print $2}')
