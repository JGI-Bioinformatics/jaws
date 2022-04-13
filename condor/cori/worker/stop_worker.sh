#!/bin/bash

echo "Stopping worker $(hostname)!"
kill $(ps aux | grep -v grep | grep -i condor_master | awk '{print $2}')
