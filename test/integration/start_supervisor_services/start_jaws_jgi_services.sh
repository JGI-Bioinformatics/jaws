#!/bin/bash

set -o pipefail

module rm python
module load python/3.8.2-dll
jaws_install="/global/home/groups-sw/lr_jgicloud/jaws-install"
gitlab_install="/global/home/groups-sw/lr_jgicloud/jaws_ci_runner"

me=$(whoami)
if [[ $me != 'jaws' ]];then
    echo "You must be running as jaws. You are $me. Please run \"collabsu jaws\"."
    exit
fi

function is_gitlab_runner_running() {
    cmd="ps -ef | grep jaws | grep gitlab-runner | grep -v grep | awk '{print \$2}'"
    pid=$(eval "$cmd")
    if [[ $pid ]]; then
        return 0
    else
        return 1
    fi
}

function start_gitlab_runner() {
    if [[ ! -e "$gitlab_install/usr/bin/gitlab-runner" ]]; then
        echo "file missing: $gitlab_install/usr/bin/gitlab-runner"
        echo "gitlab-runner failed to re-start"
        exit
    fi
    if [[ ! -e "$gitlab_install/configuration/config.toml" ]]; then
        echo "file missing: $gitlab_install/configuration/config.toml"
        echo "gitlab-runner failed to re-start"
        exit
    fi
    cmd="nohup $gitlab_install/usr/bin/gitlab-runner run -c $gitlab_install/configuration/config.toml &"

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
            echo "supervisord failed to start with command: $cmd!!!"
        fi
    fi
}

# starting services for jaws on cori

# jaws-site services
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH $branch jaws-site
    config_file="$jaws_install/jaws-supervisord-${branch}/supervisord-jaws.conf"
    start_supervisord "$config_file" "$jaws_install/jaws-supervisord-${branch}"
done

# jaws-jtm services
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH $branch jaws-jtm
    config_file="$jaws_install/jaws-supervisord-${branch}/supervisord-jtm.conf"
    start_supervisord "$config_file" "$jaws_install/jaws-supervisord-${branch}"
done

## set up monitoring services
for branch in dev staging prod; do
    echo "------------------------"
    echo ON BRANCH $branch jaws-monitor 
    config_file="$jaws_install/jaws-monitor-${branch}/supervisor/supervisord-monitor.conf"
    start_supervisord "$config_file" "$jaws_install/jaws-monitor-${branch}/supervisor"  
done

# start gitlab runner on cori
echo "------------------------"
start_gitlab_runner
