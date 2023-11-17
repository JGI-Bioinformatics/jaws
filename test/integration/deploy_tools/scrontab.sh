#!/usr/bin/env bash
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

function check_job_running {
  /usr/bin/squeue --noheader -n "jaws_${JAWS_SITE_NAME}_${JAWS_DEPLOYMENT_NAME}" --state=running,pending -u "$USER" -o "%12i %2t %9u %25j %6D %10M %12q %8f %18R"
}

function stop_worker {
  for sz in "small medium large xlarge"
  do
      /usr/bin/scancel -n "jaws_${JAWS_SITE_NAME}_${JAWS_DEPLOYMENT_NAME}_htcondor_worker_$sz" -u "$USER" 
  done
}

function stop_service {
  readarray -t job_status < <(check_job_running)
  
  # do nothing if no jobs running/pending
  if [[ ${#job_status[@]} -lt 1 ]]; then
     return
  fi

  for i in ${!job_status[@]}; do
    jobid=$(echo ${job_status[$i]} | awk '{print $1}')
    state=$(echo ${job_status[$i]} | awk '{print $2}')
    if [[ $state == "R" ]]; then
        echo "Stopping job $jobid"
        srun --overlap --jobid="$jobid" "$JAWS_BIN_DIR/supervisorctl" stop jaws-site:*
        scancel --cron "$jobid"
    else
        # if job is pending(PD) then we don't need to stop the service with srun, just cancel it.
        scancel --cron "$jobid"
    fi
  done
}

function start_service {
  if [[ ! -f "$JAWS_BIN_DIR/starter.sh" ]]; then
    echo "starter.sh not found in $JAWS_BIN_DIR, creating it now"
    envsubst '${JAWS_SITE_NAME},${JAWS_DEPLOYMENT_NAME},${JAWS_BIN_DIR},${JAWS_LOGS_DIR}' < "$DIR/templates/starter.sh" > "$JAWS_BIN_DIR/starter.sh"
  fi
  # In case our cronjob runs before we manually run our starter.sh, we attempt to kill any remaining job
  # and submit again
  stop_service
  "$JAWS_BIN_DIR/starter.sh"
}

function scrontab {
  local strategy="${1:-}"
  case "$strategy" in
    "start")
      start_service
      ;;
    "stop")
      stop_service
      stop_worker
      ;;
    *)
      echo "Unknown strategy: $strategy"
      exit 1
      ;;
  esac
}

scrontab "$@"

