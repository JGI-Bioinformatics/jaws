# Changelog for

### Mar 24, 2021

* I had MYCWD hardcoded to a path that didn't exist anymore.  It is now set by MYCWD=$(pwd).

### Mar 22, 2021

* I now check jaws status results for cori and jgi separately. If the cori services are down, but jgi is up, and I want to submit a run to jgi, then the run will continue to be submitted.

### Mar 10, 2021

* I moved the code that checks for successful jaws submission to right below the submission, since it checks the return code.
* I created another version 3.2 of the script to handle the new jaws submission command for staging (i.e. no output file needed).

### Mar 05, 2021

* I have added a better error message for when the JAWS run times out due to exhausting all the tries (i.e. CHECK_TRIES=50).  Right now the wrapper will wait for (50 x 600sec) to complete. If there are many jobs in the queue, the job may not even start running by this time and could explain this error.

### Mar 03, 2021

* I added error handling for when globus credentials are expired. A clear message is now printed to slack.

### Jan 19, 2021

* I removed alot of redundant functions.  Almost all functionality is absorbed into the wait_for_one_run function.

### Dec 07, 2020

* There was a change on how we wait for jaws runs before checking status.  Before, we just waited a specified amount of time, but when the queue is full, this time can be insufficient. So a better way was implemented that checks status and quits after success and will continue in loop for many interations, in case the queue is full.  This new function is called "wait_for_one_run".

### Dec 14, 2020

* The submit_and_wait_for_success_prod.v2.sh script was not completing 50 loops before quiting. CHECK_TRIES was set to 50, but the while loop only would go for 5 tries. For some reason, this form of the conditional (if [[ $tries < 50 ]]) was only using the first digit (i.e. 5). I had to use -lt instead of <.
