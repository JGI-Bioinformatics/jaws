#!/bin/bash -l

# jfroula slack channel (for testing)
WEBHOOK_URL=https://hooks.slack.com/services/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

# jaws
#WEBHOOK_URL=https://hooks.slack.com/services/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

# jaws_health
#WEBHOOK_URL=https://hooks.slack.com/services/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

SITE=
RELEASE=
RUNIDS_FILE="runids_${SITE}"
PROD="prod"
STAGING="staging"
CHECK_TRIES=50
CHECK_SLEEP=300

# TODO use env that is not under jfroula 
MYCWD=/global/cscratch1/sd/jaws/jfroula/test_jaws_health
JAWS_PROD_VENV=/global/cfs/projectdirs/jaws/jaws-prod/jaws-prod.sh
JAWS_STAGING_VENV=/global/cfs/projectdirs/jaws/jaws-staging/jaws-staging.sh 

usage()
{
cat<<EOF
  Usage: $0 
    <-s site [cori|jgi] (required)> 
    <-r release [prod|staging] (required)> 
EOF
exit 1
}


# check for all required args
ARGS=$*
if [[ $ARGS =~ "-s" ]] && [[ $ARGS =~ "-r" ]]; then 
    echo; 
else
    echo -e "Missing some required arguments\n"
    usage
fi


while getopts 's:r:' OPTION
do 
  case $OPTION in 
  s)    SITE="$OPTARG"
        ;;
  r)    RELEASE="$OPTARG"
        ;;
  ?)    echo usage
        exit 1
        ;;
  esac
done


function submitJob {
    RELEASE=$1
    SITE=$2
	local _runid=$3

    # submit a job that should take less than a minute once a jtm worker is scheduled.
	echo "### Submitting to $RELEASE on $SITE ###"
    jaws run submit fq_count.wdl fq_count.json ${SITE}_${RELEASE}_out_dir $SITE > ${SITE}_${RELEASE}_run_submit 
	if [[ $? > 0 ]]; then
        DATA="JAWS $RELEASE $VERSION DOWN at $SITE: jaws run submit command failed."
        echo $DATA
        curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
		exit 1
	fi

    cat ${SITE}_${RELEASE}_run_submit
    
    # capture the job id
    RUNID=$(grep '"run_id"' ${SITE}_${RELEASE}_run_submit | awk -F: '{print $2}' | tr -d ', ')
	_runid=$RUNID
}



function wait_for_one_run {
	# wait for jaws run to complete. Exit out of while loop if job completes with "succeeded" or "failed"
    RUN_ID=$1
	RELEASE=$2
    tries=1
    TMPFILE=STATUS.tmp.$$

    while [[ $tries -lt $CHECK_TRIES ]]; do
        sleep $CHECK_SLEEP

		if [[ $STATUS ]]; then
			printf "$tries try ... result is: $RESULT and status is: $STATUS\n";
		else
        	printf "$tries try...\n"
		fi

        tries=$(echo $tries + 1 | bc)

        status_result=$(jaws run status $RUN_ID > $TMPFILE)

    	if [[ $? != 0 ]]; then 
    		echo "Error \"jaws run status $RUNID > $TMPFILE\""
            DATA="JAWS $RELEASE $VERSION DOWN at $SITE: status command failed."
            echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
    		exit 1
    	fi
    
        RESULT=$(grep '"result"' $TMPFILE | awk -F: '{print $2}' | sed 's/^ *//' | tr -d ',"')
        STATUS=$(grep '"status"' $TMPFILE | awk -F: '{print $2}' | sed 's/^ *//' | tr -d ',"')
    	cromwell_id=$(grep cromwell_run_id $TMPFILE | awk -F: '{print $2}' | tr -d ' ",')
    	site=$(grep '"site_id":' $TMPFILE | awk -F: '{print $2}' | tr -d ' ",')

    	if [[ ! $RESULT ]] || [[ ! $STATUS ]] || [[ ! $cromwell_id ]] || [[ ! $site ]]; then
    		echo "Error: one or more variables could not be parsed from $TMPFILE...exiting"
    		DATA="Error: one or more variables could not be parsed from $TMPFILE...exiting"
    		echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
    		exit 1
    	fi

		# return if complete and succeeded
        if [[ $RESULT =~ "succeeded" ]] && [[ $STATUS =~ "download complete" ]]; then
            printf "Job has successfully completed\n"
			rm $TMPFILE
            DATA="JAWS $RELEASE $VERSION UP at $SITE, where run_id: [$RUN_ID] completed successfully with status [$STATUS] and result: [$RESULT]."
    		echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
            return
        fi

		# return if complete but failed
        if [[ $RESULT =~ "failed" ]] && [[ $STATUS =~ "download complete" ]]; then
            printf "Job has completed but failed. See $TMPFILE\n"
            DATA="JAWS $RELEASE $VERSION DEGRADED at $SITE. Failed to complete run, where run_id: [$RUNID]. Job completed with with status: [$STATUS] and result: [$RESULT]."
            echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
            return
        fi

    done

	# Print error if the run did not complete in the given time limit.
    printf "You have exhausted all your tries when checking jaws run status\n"

    DATA="JAWS $RELEASE $VERSION DEGRADED at $SITE. Failed to complete run, where run_id: [$RUNID]. Job completed with with status: [$STATUS] and result: [$RESULT]."
    echo $DATA
    curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
}


################ 
##### MAIN #####
################ 

# set up the jaws env
cd $MYCWD 
if [[ $RELEASE == $PROD ]]; then
	source $JAWS_PROD_VENV
elif [[ $RELEASE == $STAGING ]]; then
	source $JAWS_STAGING_VENV
else
	echo "Error: your RELEASE needs to be set to [prod|staging].  You have $RELEASE"
	exit 1
fi

VERSION=$(jaws info | awk '$1 ~ "version" {print $2}' | tr -d '"')

# test that JAWS services are up
jaws status
if [[ $? > 0 ]]; then
	DATA="JAWS services are down at $RELEASE $VERSION on $SITE"
    echo $DATA
    curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
	exit 1
fi

# submit runs to prod or staging
RUNID=
submitJob $RELEASE $SITE $RUNID  # pass an empty RUNID variable and it will get set in the function

# wait for jobs to finish & check the status
wait_for_one_run $RUNID $RELEASE

#rm -r *_run_submit *_out_dir *_runStatus
