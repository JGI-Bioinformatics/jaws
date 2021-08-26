#!/bin/bash -l

set -o pipefail

supervisor_jaws="/global/cfs/projectdirs/jaws/jaws-install"
supervisor_monitor="/global/cfs/projectdirs/jaws"
gitlab_runner="/global/cfs/cdirs/m342/jaws_runner/usr/bin/gitlab-runner"

myhost=$(hostname)
me=$(whoami)

if [[ $myhost != 'cori20' ]];then
    echo "Host must be cori20. You are on $myhost. Please run \"ssh cori20\"."
    exit
fi
if [[ $me != 'jaws' ]];then
    echo "You must be running as jaws. You are $me. Please run \"collabsu jaws\"."
    exit
fi

## Function to get supervisord pid
function is_gitlab_runner_running() {
    cmd="ps -ef | grep jaws | grep gitlab-runner | grep -v grep | awk '{print \$2}'"
    echo "$cmd"
    pid=$(eval "$cmd")
    if [[ $pid ]]; then
        return 0
    else
        return 1
    fi

}

function start_gitlab_runner() {
    if [[ ! -e $gitlab_runner ]]; then
        echo file missing or empty: $gitlab_runner
        echo gitlab-runner failed to re-start
        exit
    fi
    cmd="nohup ./$gitlab_runner run &"
    if is_gitlab_runner_running; then
        echo gitlab runner is already running
    else
        echo re-starting gitlab runner 
        echo "$cmd"
        eval "$cmd"

        if is_gitlab_runner_running; then
            echo gitlab runner started successfully
        else
            echo gitlab runner failed to restart!!!!
        fi

    fi
}


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

# starting services for jaws on cori
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH $branch jaws-site
    config_file="$supervisor_jaws/jaws-supervisord-${branch}/supervisord-jaws.conf"
    start_supervisord "$config_file" "$supervisor_jaws/jaws-supervisord-${branch}"
done

## set up monitoring services
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH $branch jaws-monitoring
    config_file="$supervisor_monitor/jaws-monitor-${branch}/supervisor/supervisord-monitor.conf"
    start_supervisord "$config_file" "$supervisor_monitor/jaws-monitor-${branch}/supervisor"  
done

# start gitlab runner on cori 
echo "------------------------"
start_gitlab_runner
