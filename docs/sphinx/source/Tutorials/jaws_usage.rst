======================== 
How to Run a WDL in JAWS
========================

.. role:: bash(code)
   :language: bash

.. warning::
   Before these commands will work, you need to set up everything:
   See `Quickstart Example <jaws_quickstart.html>`_ for full setup instructions.


Example submitting a workflow:

:bash:`jaws run submit <wdl> <inputs json> <out dir> <site>` 

... where "run" is a command, and "submit" is a sub-command. Running each command/sub-command without arguments 
shows the help message.


*************
JAWS commands
*************


There are four top level JAWS commands (but many sub commands):

  :bash:`jaws login` do once to obtain a Globus token. Just follow directions after running the command. 

  :bash:`jaws run` command for job management (e.g. submission, monitoring and log-viewing). 

  :bash:`jaws wdl`  command for workflow management (e.g. add, delete & view existing workflows). 

  :bash:`jaws status`  displays the status of jaws-services at different sites. 


*jaws login*:
----------------------

You only need to run this once to get a Globus token.  Follow the directions to obtain the token from Globus and this will allow JAWS to transfer data around via Globus in your name.

.. code-block:: bash

    jaws login
    
    # This will generate a URL that you need to follow to get the Globus token. Then save the token according to directions. 
    


*jaws run*:
-------------------

.. code-block:: text

  cancel       cancel a run; prints whether aborting was successful or not.
  history      print a list of the user's past runs, including failed runs.
  inputs       generate inputs template (JSON) from workflow (WDL) file.
  list-sites   list available Sites to run your jobs
  log          view the log of Run state transitions. this is at a higher level than task-log.
  metadata     print the detailed metadata for a run, partly produced by Cromwell server.
  output       view the stdout/stderr output of Tasks.
  queue        list user's unfinished runs.
  status       print the current status of a run.
  submit       submit a run for execution at a JAWS-Site.
  task-log     get log of each Tasks' state transitions.
  task-status  show the current status of each task.
  validate     validate a WDL using Cromwell's WOMTool.

   
*jaws wdl*:
-------------------

.. code-block:: text

  about       return README document for a workflow
  add         add a workflow to the catalog. do this before **release**
  delete      remove a workflow from the catalog
  get         get WDL specification for a workflow. this generates a WDL for running.
  list        list shared workflows in the catalog
  release     mark a version as an immutable production release. the README can still be changed.
  update-doc  update a workflow README in the catalog
  update-wdl  update a workflow WDL in the catalog (unavailable for released WDLs)
  versions    list available versions of a workflow for the public to run. some versions may not be publicly available.
 

*jaws status*:
----------------------

This command shows the status of the various jaws-services. Some services are site specific.

.. code-block:: text

    {
      "JAWS-Central": "UP",
      "JGI-Cromwell": "Unknown",
      "JGI-RMQ": "UP",
      "JGI-Site": "DOWN",
      "CORI-Cromwell": "UP",
      "CORI-RMQ": "UP",
      "CORI-Site": "UP"
    }



Examples
--------

**To run a wdl**

.. code-block:: text

    # find available sites and submit to CORI
    jaws run list-sites
    jaws run submit my.wdl inputs.json out cori

    # submit it to JGI
    jaws run submit my.wdl inputs.json out jgi


**Anyone can share a WDL. To see a list of workflows available in the catalog run**

.. code-block:: text

  jaws wdl list

  # output: where fq_count is the name of the WDL and dev is the version.  
  [
      "fq_count",
      "dev",
      "ekirton",
      "2020-03-24T02:04:10Z",
      "2020-03-24T09:14:18Z"
  ]


**To see info about that workflow (generated from a README)**

.. code-block:: text

   # note that a version is required
   jaws wdl about fq_count dev 


**To run a WDL from the catalog, there are a couple extra steps (from "jaws run list" we saw there is a wdl in the catalog called fq_count)**

.. code-block:: text

    # create the wdl
    jaws wdl get fq_count dev > my.wdl
    
    # create a template for inputs.json 
    jaws run inputs my.wdl > inputs.json

    # cusomize the values in inputs.json
    vi inputs.json

    # run as usual
    jaws run submit my.wdl inputs.json out cori


.. note::

    From any job submission, you can see a run id (i.e. see 121 below). You'll use this for future commands.

.. code-block:: text

  # output looks like
  {
  "output_dir": "<full_path>/out",
  "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
  "run_id": 121,
  "site_id": "CORI",
  "status": "uploading",
  "submission_id": "7d2606b9-569f-4d50-9423-c1acb5441c6b",
  "upload_task_id": "07ffa460-88ac-11ea-b3ba-0ae144191ee3"
  }



**See the status of a Run using job ID**

.. code-block:: text

  jaws run status 121


**Monitoring Runs**

When monitoring the runs, each task transitions between the following states. 

.. code-block:: text

   uploading            # input data are being copied to scratch by Globus
   missing input        # run was uploaded but some of the required files were missing
   upload complete      # Globus finished copying all your files to scratch
   submitted            # job submitted to JTM and worker pools have been requested
   queued               # waiting for worker pools to be reserved from cluster
   running              # the run is being executed by Cromwell
   succeeded            # Cromwell completed the run but results need to be transfered
   ready                # results are ready for Globus transfer off of site scratch
   downloading          # results are being copied by Globus
   download complete    # results have been copied to your output directory. signifies end of run
   failed               # runing error from either jaws or user's wdl
   canceled             # run was cancelled by user or JTM issue


**Checking the staus of a task**

.. code-block:: text

    # the two status commands show the current status of the run or tasks of the run
    jaws run status 121
    jaws run task-status 121

    # the log commands show all the past states of either the run or tasks of the run
    jaws run log 121
    jaws run task-log 121

**Get current or old history of jobs owned by you**

.. code-block:: text

   # get list of your currently running jobs
   jaws run queue                                      
   
   # view history of your jobs for last 7 days 
   jaws run history --days 7


