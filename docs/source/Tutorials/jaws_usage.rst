****************
How to Run a WDL
****************
Three ways to run a WDL:

.. role:: bash(code)
   :language: bash

1) using the :bash:`cromwell.jar` command - when testing and developing a WDL.  
2) using the :bash:`jaws` command - when you are ready, you can run the WDL file within JAWS.  
3) using the :bash:`jaws` command - if there is a WDL registered in JAWS (i.e. by someone else), you can run it.


Running Outside JAWS
--------------------
Do this when you are developing a WDL. You will be running all your tasks locally but turn-around should be faster. 

.. code-block:: bash

   ssh cori
   java -jar /global/dna/projectdirs/DSI/workflows/cromwell/java/cromwell.jar run <your.wdl> --inputs <input.json>


You can run your tasks as sbatch (slurm) jobs but you'll have to run with a config file.  See :doc:`Configuring Backends </Specifications/configure_cromwell>`

Running Inside JAWS
-------------------
There are two ways to run a WDL

   1) use a WDL file from your **working directory** 

   2) use a WDL that is **registered** in the JAWS catalog. 

The command is the same in both cases :bash:`jaws submit <wdl> <inputs>`

To see existing workflows in the catalog, run :bash:`jaws list`


*************
JAWS commands
*************

.. note::
   Before these commands will work, you need to install the conda environment:
   "conda activate jaws". Please see :doc:`Set up the JAWS environment </Specifications/setting_up_environment>`

There are two main JAWS commands:

  :bash:`jaws` command deals with jobs (e.g. submission, monitoring and log-viewing) 

  :bash:`wf`  command deals with workflows (e.g. description of & listing of) 

optional arguments:
-------------------

**jaws command**

.. code-block:: bash

   JAWS : JGI Analysis Workflows Service
   -------------------------------------
   
   jaws submit <wf file or name> <inputs> (subs)    : submit a job using a workflow name or WDL file
   jaws batch <wf file or name> <inputs> (subs)     : submit a batch of jobs
   
   jaws cancel <job_id>                             : cancel a job
   
   jaws queue                                       : list your unfinished jobs
   jaws last <num_days> (wf_name)                   : list all recent jobs, optionally limited to a workflow
   
   jaws status <job_id>                             : current status of a job
   jaws tasks <job_id>                              : current status of each task of a job
   
   jaws log <job_id>                                : brief information about a job
   jaws metadata <job_id>                           : detailed information about a job
   jaws outputs <job_id>                            : list output files of a successfully completed job
   
   jaws invalidate <job_id> (task_name)             : prevent a job outut from being reused by cache
   jaws wait <job_id>                               : block until job is done (rc indicates success or failure). This is useful when calling
													: jaws within a bash script.
   

**wf command**

.. code-block:: bash

   JAWS Workflows Catalog
   ----------------------

   wf list                                        : list shared workflows
   wf about <name/version>                        : display readme doc for a workflow
   wf wdl <name/version>                          : return the WDL for a workflow
   wf inputs <name/version>                       : generate inputs template


Examples
--------

to see a list of workflows

::

  wf list

  # output: where bbstats is the name of the WDL and 1.0.0 is the version.  
  {
    "bbstats/1.0.0": "http://app.jaws-svc.prod-cattle.stable.spin.nersc.org:60045/api/workflows/bbstats",
    ...
  }



to see info about that workflow

::

   # note that no version is required here
   wf about bbstats

to create a template for your inputs file (e.g. inputs.json).

::

   wf inputs bbstats/latest


to submit a job 

::

  # use registered wdl from the above list (you need to supply the inputs.json; 
  # or test with /global/project/projectdirs/jaws/jgi-workflows/bbstats/test.json)
  jaws bbstats/latest inputs.json
 
  # supply your own wdl and inputs
  jaws submit 1.1.0.wdl inputs.json


or to see the status or metadata of a run using job ID

::

  job status ec43alkoi22342kloiaudkjo909ad

  # there's alot of good stuff in metadata so check it out
  job metadata ec43alkoi22342kloiaudkjo909ad


See log from cromwell

::

   job log ec43alkoi22342kloiaudkjo909ad


get current or old history of jobs

::

   # get list of your currently running jobs
   jaws queue                                      
   
   # view history of your jobs for last 7 days 
   jaws last 7 
   jaws last 7 bbstats/1.0.0


clear cache

Use this when you want to re-run one or more of your tasks in your workflow (i.e. don't use cached results).
For example, if you change something in a script but the WDL doesn't change, you will use cached results (which will not reflect changes in your script).

::

   jaws invalidate ec43alkoi22342kloiaudkjo909ad trimFastq
   # now re-submit the wdl to jaws.
   jaws submit 1.1.0.wdl inputs.json


Additional Commands
-------------------

**UTILITIES**

:bash:`wfcopy.py <src_dir> <dest_dir>`  ---  copy and flatten a job's output to another dir
