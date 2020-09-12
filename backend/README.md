## JAWS Backend

This package provides a simple command-line application to enable Cromwell communication with jaws-site via RPC.

## Cromwell Backend

Cromwell's flexible configuration allows one to use any compute resource that can be accessed via the
command-line, described below.  Common backends include cluster schedulers (e.g. SLURM, SGE, Condor) as well as
cloud computing resources (e.g. Google Cloud, AWS).  Here we provide our own middleware between Cromwell and
cluster schedulers in order to add several key features.

## Operations

Cromwell requires only three operations of a backend:

  * submit : submit a task for execution given a shell script; returns job_id
  * kill : abort execution of the script, indicated by job_id
  * check_status : check if a task is running, indicated by job_id

No other operations are required or supported.

## Features

The following features are made possible by inserting the JAWS middleware between Cromwell and the scheduler:

  * priority queue : process tasks from earlier runs first, regardless of order of submission
  * node reuse : process multiple tasks per job reservation, reducing queue-wait and load on scheduler
