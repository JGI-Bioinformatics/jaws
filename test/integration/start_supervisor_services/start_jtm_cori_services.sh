#!/bin/bash

set -o pipefail

supervisor_install="/tmp"

myhost=$(hostname)
me=$(whoami)

if [[ $myhost != 'cori20' ]];then
    echo "Host must be cori20. You are on $myhost. Please run \"ssh cori20\"."
    exit
fi
if [[ $me != 'jaws_jtm' ]];then
    echo "You must be running as jaws_jtm. You are $me. Please run \"collabsu jaws_jtm\"."
    exit
fi

## Function to get supervisord pid
function is_supervisord_running() {
    local config_file=$1
    local supervisor_dir=$2
    cmd="$supervisor_dir/bin/supervisorctl -c $config_file status"
    eval "$cmd" 
    if [[ $? == 4 ]]; then
        return 1
    else
        return 0
    fi
}


## Function to start supervisord
function start_supervisord() {
    local config_file=$1
    local supervisor_dir=$2
    if is_supervisord_running "$config_file" "$supervisor_dir"; then
        echo supervisord already running
    else
        echo "Starting supervisord ..."
        cmd="$supervisor_dir/bin/supervisord -c $config_file"
        echo "$cmd"
        eval "$cmd"

        if is_supervisord_running "$config_file" "$supervisor_dir"; then
            echo "supervisord started successfully"
        else
            echo "supervisord failed to start with command: $cmd"
        fi
    fi
}

# starting services for jaws_jtm on cori
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH $branch jaws-jtm
    config_file="$supervisor_install/jaws-supervisord-${branch}/supervisord-jtm.conf"
    start_supervisord "$config_file" "$supervisor_install/jaws-supervisord-${branch}"
done
