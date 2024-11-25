#!/bin/bash

BUFFER_MIN=120 # Number of minutes to buffer before maintenance
SCRONTAB_NAME="jaws_${JAWS_SITE_NAME}_${JAWS_DEPLOYMENT_NAME}" # Name of the scrontab job to check

MY_JOB=$(/usr/bin/squeue --noheader -n "$SCRONTAB_NAME" --state=running,pending -u "$USER" -o "%12i %2t %9u %25j %6D %10M %12q %8f %18R")

if [[ "${#MY_JOB}" -ge 2 ]]; then
    printf "[%s] Job %s already in the queue \n" "$(date "+%m-%d-%Y-%H:%M:%S")" "$SCRONTAB_NAME"
    printf "===============================\n%s\n===============================\n" "$MY_JOB"
    exit 0;
fi

STR_TIME=$(scontrol show res -o | egrep 'login|workflow' | awk -F'[= ]' '{print $4}'| head -n 1)

mkdir -p ${JAWS_LOGS_DIR}/site-slurm
sbatch_cmd="sbatch --job-name=${SCRONTAB_NAME} --dependency=singleton --time=90-00:00:00 --output=${JAWS_LOGS_DIR}/site-slurm/jaws-site-slurm-%j.out --error=${JAWS_LOGS_DIR}/site-slurm/jaws-site-slurm-%j.err --cpus-per-task=16"

if [[ "$STR_TIME" == "" ]]; then
    printf '[%s] No maintenance found\n' "$(date "+%m-%d-%Y-%H:%M:%S")"
else
    BAD_TIME=$(date "+%s" -d "$STR_TIME")
    CUR_TIME=$(date "+%s")
    TIME_LEFT_MIN=$(((BAD_TIME-CUR_TIME)/60))

    printf '[%s] Minutes before next maintenance: %s\n' "$(date "+%m-%d-%Y-%H:%M:%S")" "$TIME_LEFT_MIN"

    TIME_MIN=$((TIME_LEFT_MIN-BUFFER_MIN))
    if [[ "$TIME_MIN" -gt 0 ]]; then
      sbatch_cmd+=" --time-min=${TIME_MIN}"
    fi
fi


sbatch_cmd+=" --qos workflow -C cron ${JAWS_BIN_DIR}/jaws-site-cronjob"
unset ${!SLURM_@};
echo "[$(date "+%m-%d-%Y-%H:%M:%S")] $sbatch_cmd"
eval "$sbatch_cmd"
