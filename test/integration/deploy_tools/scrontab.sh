#!/usr/bin/env bash

function stop_service {
  squeue --format="%.18i %.9P %.35j %.8u %.8T %.10M %.9l %.6D %R" --me | grep "${JAWS_SITE_NAME}_${JAWS_DEPLOYMENT_NAME}" |  xargs -n 1 scancel || true
}

function start_service {
  echo "Login to perlmutter and modify the scrontab with scrontab -e"
}

function scrontab {
  local strategy="${1:-}"
  case "$strategy" in
    "start")
      start_service
      ;;
    "stop")
      stop_service
      ;;
    *)
      echo "Unknown strategy: $strategy"
      exit 1
      ;;
  esac
}

scrontab "$@"


