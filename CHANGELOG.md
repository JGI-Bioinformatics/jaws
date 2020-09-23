# Changelog

## 2.1

 - We no longer reformat output automatically; instead users may use "wfcopy" command to reformat the output after their run completes.  Users may find it useful for runs performed by Cromwell, outside of JAWS.
 - The output of tasks are returned as the tasks complete, rather than waiting for the end.  As the output includes the stdout/stderr files, the "run outputs" command has been deprecated and an new "run errors" command summarizes any errors (replaces "run output --failed").
 - "run task-log" and "run task-status" are now real-time; previously were out of sync by up to 10 seconds due to update interval.
 - Simplified user config file, provided a .sh file to source for activating jaws which simplifies use of multiple jaws deployments, and provide wheel file if you wish to install your own client or use it in your own python software.
 - Add "info" command to provide jaws version and deployment user client is using, as well as the link to the documentation appropriate for that release.
 - "run list-sites" now includes the maximum requestable RAM available at each Site
 - "run metadata" now returns any subworkflows' metadata too.
 - assorted minor bugfixes and improvements to usability.