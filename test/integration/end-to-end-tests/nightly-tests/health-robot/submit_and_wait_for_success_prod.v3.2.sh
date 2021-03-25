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
CHECK_SLEEP=600

# TODO use env that is not under jfroula 
MYCWD=$(pwd)
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

# make sure site is all lower case
SITE=$(echo "$SITE" | awk '{print tolower($0)}')

function submitJob {
    RELEASE=$1
    SITE=$2
	local _runid=$3

    # submit a job that should take less than a minute once a jtm worker is scheduled.
	echo "### Submitting to $RELEASE on $SITE ###"
    jaws run submit fq_count.wdl fq_count.json $SITE > ${SITE}_${RELEASE}_run_submit 2> ${SITE}_${RELEASE}.stderr

	if [[ $? > 0 ]]; then
		ERR_MSG=$(cat ${SITE}_${RELEASE}.stderr)
        DATA="JAWS $RELEASE $VERSION submit command failed. $ERR_MSG"
        echo $DATA
        curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
		exit 1
	fi

	# Test that my globus credentials are active and return an appropriate error if they are not.
	grep_results=$(grep 'Globus endpoint has expired' ${SITE}_${RELEASE}.stderr)
	if [[ $grep_results ]]; then
		DATA="JAWS ${RELEASE} on ${SITE}: (Globus credentials need to be renewed!)"
		cat ${SITE}_${RELEASE}.stderr
        curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
		exit 1
	fi

    cat ${SITE}_${RELEASE}_run_submit
    
    # capture the job id
    RUNID=$(grep '"run_id"' ${SITE}_${RELEASE}_run_submit | awk -F: '{print $2}' | tr -d ', ')
if [[ $? > 0 ]]; then
        DATA="JAWS $RELEASE $VERSION failed to parse the RUNID from file: ${SITE}_${RELEASE}_run_submit ...exiting"
        echo $DATA
        curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
		exit 1
	fi

	rm ${SITE}_${RELEASE}.stderr ${SITE}_${RELEASE}_run_submit

	# set _runid to set a global variable
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

        status_result=$(jaws run status $RUN_ID > $TMPFILE 2> ${RUN_ID}_status_out)

    	if [[ $? != 0 ]]; then 
			ERR_MSG=$(cat ${RUN_ID}_status_out)
            DATA="JAWS release: $RELEASE $VERSION DOWN at $SITE: status command failed.  $ERR_MSG"
            echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
    		exit 1
    	fi
    
        RESULT=$(grep '"result"' $TMPFILE | awk -F: '{print $2}' | sed 's/^ *//' | tr -d ',"')
        STATUS=$(grep '"status"' $TMPFILE | awk -F: '{print $2}' | sed 's/^ *//' | tr -d ',"')
    	cromwell_id=$(grep cromwell_run_id $TMPFILE | awk -F: '{print $2}' | tr -d ' ",')
    	site=$(grep '"site_id":' $TMPFILE | awk -F: '{print $2}' | tr -d ' ",')

    	if [[ ! $RESULT ]] || [[ ! $STATUS ]] || [[ ! $cromwell_id ]] || [[ ! $site ]]; then
    		DATA="Error: one or more variables could not be parsed from $TMPFILE...exiting"
    		echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
    		exit 1
    	fi

		# return if complete and succeeded
        if [[ $RESULT =~ "succeeded" ]] && [[ $STATUS =~ "download complete" ]]; then
			rm $TMPFILE
            printf "Job has successfully completed\n"
            DATA="JAWS $RELEASE $VERSION UP at $SITE, where run_id: [$RUN_ID] completed successfully with status [$STATUS] and result: [$RESULT]."
    		echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
            return
        fi

		# return if complete but failed
        if [[ $RESULT =~ "failed" ]] && [[ $STATUS =~ "download complete" ]]; then
			rm $TMPFILE
            printf "Job has completed but failed. See $TMPFILE\n"
            DATA="JAWS $RELEASE $VERSION DEGRADED at $SITE. Failed to complete run, where run_id: [$RUNID]. Job completed with with status: [$STATUS] and result: [$RESULT]."
            echo $DATA
            curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
            return
        fi

    done

	# Print error if the run did not complete in the given time limit.
	rm $TMPFILE
    printf "You have exhausted all your tries when checking jaws run status\n"

    DATA="JAWS $RELEASE $VERSION TIMED OUT at $SITE. Failed to complete run because number of tries [$CHECK_TRIES] have been exhausted., where run_id: [$RUNID]. Job completed with with status: [$STATUS] and result: [$RESULT]. This error could be caused by long wait times in the queue."
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
# Need to check both cori and jgi services separately. For example,
# just because cori services are down doesn't mean jgi is down too.
#
# "CORI-Cromwell": "UP",
# "CORI-RMQ": "UP",
# "CORI-Site": "UP",
# "JAWS-Central": "UP",
# "JGI-Cromwell": "UP",
# "JGI-RMQ": "UP",
# "JGI-Site": "UP"

# test cori services
TMPFILE_SERVICES=services.tmp.$$
jaws status > $TMPFILE_SERVICES 2>&1

CORI_SERVICES_RESULTS=$(grep DOWN $TMPFILE_SERVICES | grep CORI)
if [[ $CORI_SERVICES_RESULTS ]]; then

    # we only need to report error and quit if we're trying to run at cori.
    if [[ $SITE == "cori" ]]; then
        DATA="JAWS services are down for CORI $RELEASE $VERSION"
        echo $DATA
        curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
        exit 1
    else
        echo "we were not trying to use cori so I'll continue submitting run."
    fi
else
    echo "cori services are all up."
fi

# test jgi services
JGI_SERVICES_RESULTS=$(grep DOWN $TMPFILE_SERVICES | grep JGI)
if [[ $JGI_SERVICES_RESULTS ]]; then

    # we only need to report error and quit if we're trying to run at jgi.
    if [[ $SITE == "jgi" ]]; then
        DATA="JAWS services are down for JGI $RELEASE $VERSION"
        echo $DATA
        curl -X POST -H 'Content-type: application/json' --data '{"text":"'"$DATA"'"}' $WEBHOOK_URL
        exit 1
    else
        echo "we were not trying to use JGI so I'll continue submitting run."
    fi
else
    echo "JGI services are all up."
fi

rm $TMPFILE_SERVICES

# submit runs to prod or staging
RUNID=
submitJob $RELEASE $SITE $RUNID  # pass an empty RUNID variable and it will get set in the function

# wait for jobs to finish & check the status
if [[ $RUNID ]]; then
	wait_for_one_run $RUNID $RELEASE
else
	echo "No runid found for jaws submission....exiting"
	exit 1
fi

