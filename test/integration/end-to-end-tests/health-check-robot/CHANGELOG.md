# Changelog for

### Jan 19, 2021
Version 3.0.0

* I removed alot of redundant functions.  Almost all functionality is absorbed into the wait_for_one_run function.

### Dec 07, 2020
Version 2.0.0

* There was a change on how we wait for jaws runs before checking status.  Before, we just waited a specified amount of time, but when the queue is full, this time can be insufficient. So a better way was implemented that checks status and quits after success and will continue in loop for many interations, in case the queue is full.  This new function is called "wait_for_one_run".

### Dec 14, 2020

* The submit_and_wait_for_success_prod.v2.sh script was not completing 50 loops before quiting. CHECK_TRIES was set to 50, but the while loop only would go for 5 tries. For some reason, this form of the conditional (if [[ $tries < 50 ]]) was only using the first digit (i.e. 5). I had to use -lt instead of <.
