# Summary
This repo contains scripts and supporting data for testing the health of JAWS. The script(s) in this repository are run by a cronjob on cori20 under user jfroula.

For now, the script runs a simple WDL that just counts the number of fastq reads "fq_counter" which is in the JAWS catalog.  The test is whether or not the WDL completes in 10 minutes or not, which should be enough time but depends on the slurm schedular.  Better and more tests will be added to this repo as we move ahead.

## crontjob
The submit_and_wait_for_success.sh script tests jaws-dev and runs nightly at 4:00am and sends email to jlfroula@lbl.gov as well as a notification to #jaws slack channel.
The submit_and_wait_for_success_prod.sh script tests jaws-prod and runs nightly at 4:00am and sends email to jlfroula@lbl.gov as well as a notification to #jaws slack channel.

## how to notify slack channel
The magic is this command which sends a message to slack #jaws channel.  
```
curl -X POST -H 'Content-type: application/json' --data '{"text":"JAWS failed to complete job '$JOBID'.  Job hung in state: '$STATUS'"}' $WEBHOOK_URL
```

This curl command was set up by going to 
`https://api.slack.com/apps` and creating an app, configured it and then by clicking on the "Incoming Webhooks" tab, I could add a webhook which gave me the above curl command.

## What is an end-to-end test
In the widely-used testing pyramid strategy, end-to-end tests act more like a safeguard: most of your code should be covered by unit tests that allow you to easily identify the source of a problem, should one occur. Rather, you will likely want to limit the number of end-to-end tests to just enough to 
* give you the confidence that the deployment went as intended
* that your infrastructure is up and running 
* that your units of code work well together.
