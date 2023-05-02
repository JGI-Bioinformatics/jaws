#!/usr/bin/env bash

SERVICES=("rpc-server" "run-daemon" "transfer-daemon" "perf-metrics-daemon" "task-log");

function disable_service {
  systemctl --user disable "$1"
}

function enable_service {
  systemctl --user enable "$1"
}

function start_service {
  systemctl --user start "$1"
}

function stop_service {
  systemctl --user stop "$1"
}

function stop_services {
  for service in "${SERVICES[@]}"; do
    deploy_name="$service-${JAWS_DEPLOYMENT_NAME}"
    local status=$(systemctl --user is-active "$deploy_name")
    if [[ "$status" == "inactive" ]]; then
      echo "Service $deploy_name is already stopped"
    else
       stop_service "$deploy_name"
    fi
  done
}

function start_services {
  daemon-reload

  for service in "${SERVICES[@]}"; do
    deploy_name="$service-${JAWS_DEPLOYMENT_NAME}"
    status=$(systemctl --user is-enabled "$deploy_name")
    if [[ "$status" == "disabled" ]]; then
      enable_service "$deploy_name"
    fi
    start_service "$deploy_name"
  done
}

function enable_services {
  for service in "${SERVICES[@]}"; do
    deploy_name="$service-${JAWS_DEPLOYMENT_NAME}"
    enable_service "$deploy_name"
  done
}

function disable_services {
  for service in "${SERVICES[@]}"; do
    deploy_name="$service-${JAWS_DEPLOYMENT_NAME}"
    disable_service "$deploy_name"
  done
}

function daemon-reload {
  systemctl --user daemon-reload
}

function status {
  for service in "${SERVICES[@]}"; do
    deploy_name="$service-${JAWS_DEPLOYMENT_NAME}"
    systemctl --user status "$deploy_name"
  done
}

function systemd {
  strategy="${1:-}"
  case "$strategy" in
    "start")
      start_services
      ;;
    "stop")
      stop_services
      ;;
    "enable")
      enable_services
      ;;
    "disable")
      disable_services
      ;;
    "status")
      status
      ;;
    *)
      echo "Unknown strategy: $strategy"
      exit 1
      ;;
  esac
}

systemd "$@"
