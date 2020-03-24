========================
How to Run a WDL in JAWS
========================

.. role:: bash(code)
   :language: bash

.. note::
   Before these commands will work, you need to install the conda environment:

   on CORI
     :bash:`source activate /global/cfs/projectdirs/jaws/prod/cli/`

   on LAB IT
     :bash:`source activate /global/home/groups-sw/lr_jgicloud/dev/jaws_client/`


Then you can submit a run

:bash:`jaws run submit <wdl> <inputs.json>` 


*************
JAWS commands
*************


There are three flavors of JAWS commands:

  :bash:`jaws run` command for job management (e.g. submission, monitoring and log-viewing) 

  :bash:`jaws wdl`  command for workflow catalog (e.g. view existing workflows) 

  :bash:`jaws util`  workflow developer utilities

jaws run options:
-------------------

.. code-block:: bash

   block <run_id>                 block until run is complete. 
   cancel <run_id>                cancel a run
   convert <run_id>               covert table of inputs to JSON format for submitting a batch of...
   delete <run_id>                delete the output of a run or task to avoid caching. Use when you don't want
								  to re-use any outputs from the last run.
   history <num_days> (wf_name)   list past runs for a given number of days
   metadata <run_id>              detailed information about a run
   queue                          list your unfinished runs
   status <run_id>                show current status of a run
   submit <wdl> <inputs>          submit a run for execution
   tasks <run_id>                 show status of each task of a run
   wait <run_id>                  wait until run is complete; check return code. I.e. In a script, you can run 
								  a jaws job as if it were just another bash command. You would use **wait** 
								  to let jaws complete before the next command

   

jaws wdl options:
-------------------

.. code-block:: bash

  about       return README document for a workflow
  add         add a workflow to the catalog. do this before **release**
  release     mark a version as immutable production release
  delete      remove a workflow from the catalog
  get         get WDL specification for a workflow
  list        list shared workflows
  update-doc  update a workflow README in the catalog
  update-wdl  update a workflow WDL in the catalog
  versions    list available versions of a workflow
 

jaws util options:
-------------------

.. code-block:: bash

  inputs    generate inputs template from WDL file
  status    current system status
  validate  validate your WDL
  wfcopy    copy cromwell output to specified dir


Examples
--------

to see a list of workflows

::

  jaws wdl list

  # output: where bbstats is the name of the WDL and 1.0.0 is the version.  
  {
    "bbstats/1.0.0": "http://app.jaws-svc.prod-cattle.stable.spin.nersc.org:60045/api/workflows/bbstats",
    ...
  }



to see info about that workflow

::

   # note that no version is required here
   jaws wdl about bbstats

to create a template for your inputs file (e.g. inputs.json).

::

   jaws wdl inputs bbstats/latest


to submit a job 

::

  # use registered wdl from the above list (you need to supply the inputs.json; 
  # or test with /global/project/projectdirs/jaws/jgi-workflows/bbstats/test.json)
  jaws run submit metagenome_assembly.wdl inputs.json
 

or to see the status or metadata of a run using job ID

::

  jaws run status ec43alkoi22342kloiaudkjo909ad

  # there's alot of good stuff in metadata so check it out
  jaws run metadata ec43alkoi22342kloiaudkjo909ad


get current or old history of jobs

::

   # get list of your currently running jobs
   jaws run queue                                      
   
   # view history of your jobs for last 7 days 
   jaws run history 


clear cache

Use this when you want to re-run one or more of your tasks in your workflow (i.e. don't use cached results).
For example, if you change something in a script but the WDL doesn't change, you will use cached results (which will not reflect changes in your script).

::

   jaws run delete ec43alkoi22342kloiaudkjo909ad

   # now re-submit the wdl to jaws.
   jaws run submit metagenome_assembly.wdl inputs.json


This gives us a template for "inputs.json". 

::

	jaws util inputs metagenome_assembly.wdl > inputs.json
