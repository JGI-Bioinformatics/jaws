#!/bin/bash

set -o pipefail

supervisor_dir="/opt/jaws"

myhost=$(hostname)
me=$(whoami)

if [[ $myhost != 'jaws' ]];then
    echo "Host must be jaws You are on $myhost. Please run \"ssh -l <username>@lbl.gov jaws.lbl.gov\" and use your AD password (no MFA)."
    exit
fi
if [[ $me != 'jaws' ]];then
    echo "You must be running as jaws. You are $me. Please run \"sudo -i -u jaws\"."
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

# starting services for jaws_jtm on central
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH "$branch" jaws-central
    config_file="$supervisor_dir/jaws-central-supervisord-${branch}/supervisord-jaws-central.conf"
    start_supervisord "$config_file" "$supervisor_dir/jaws-central-supervisord-${branch}"
done

# starting jaws-monitoring service on central
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH "$branch" jaws-monitor
    config_file="$supervisor_dir/jaws-monitor-${branch}/supervisor/supervisord-monitor.conf"
    start_supervisord "$config_file" "$supervisor_dir/jaws-monitor-${branch}/supervisor"
done
