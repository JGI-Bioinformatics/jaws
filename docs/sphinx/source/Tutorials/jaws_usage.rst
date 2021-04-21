======================== 
How to Run a WDL in JAWS
========================

.. role:: bash(code)
   :language: bash

.. warning::
   Before these commands will work, you need to set up everything:
   See `Quickstart Example <jaws_quickstart.html>`_ for full setup instructions.


Example submitting a workflow:

:bash:`jaws run submit --tag <give this workflow a name> <wdl> <inputs json> <site>` 

... where "run" is a command, and "submit" is a sub-command. Running any command/sub-command without arguments 
shows the help message.

.. note:: 
	using :bash:`--tag` is useful to keep track of things when you have multiple runs. The tag identifier you give will show up with the various jaws logging and status commands.


*************
JAWS commands
*************


There are four top level JAWS commands (but many sub commands):

  :bash:`jaws login` do once to obtain a Globus token. Just follow directions after running the command. 

  :bash:`jaws run` command for job management (e.g. submission, monitoring and log-viewing). 

  :bash:`jaws wdl`  command for workflow management (e.g. add, delete & view existing workflows). 

  :bash:`jaws status`  displays the status of jaws-services at different sites. 


The command: *jaws login*
-----------------------------

You only need to run this once to get a Globus token.  Follow the directions to obtain the token from Globus and this will allow JAWS to transfer data around via Globus in your name.

.. code-block:: text

    jaws login
    
    # This will generate a URL that you need to follow to get the Globus token. Then save the token according to directions. 
    


The sub-commands for *jaws run*:
--------------------------------

.. code-block:: text

  cancel       Cancel a run; prints whether aborting was successful or not.
  errors       View error messages and stderr for failed tasks.
  get          Copy the output of a run to the specified folder.
  history      Print a list of the user's past runs.
  inputs       Generate inputs template (JSON) from workflow (WDL) file.
  list-sites   List available Sites
  log          View the log of Run state transitions. this is at a higher level than task-log.
  metadata     Print the detailed metadata for a run.
  queue        List user's unfinished runs.
  status       Print the current status of a run.
  submit       Submit a run for execution at a JAWS-Site.
  task-log     Get log of each Tasks' state transitions.
  task-status  Show the current status of each task.
  validate     Validate a WDL using Cromwell's WOMTool.

The sub-commands for *jaws wdl*:
--------------------------------

.. code-block:: text

  about       Return README document for a workflow.
  add         Add a new workflow to the catalog.
  delete      Remove a workflow from the catalog.
  get         Get workflow specification (WDL) for a workflow. This generates a WDL for running.
  list        List available workflows in the Catalog.
  release     Mark a version as released, which makes it's WDL immutable. The README can still be changed.
  update-doc  Update a workflow's README in the catalog.
  update-wdl  Update a workflow's WDL in the catalog (unavailable for released WDLs)
  versions    List publicly available versions of a specified workflow. Not all versions are publicly available.

The command: *jaws status*
--------------------------

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

    # submit it to JGI
    jaws run submit --tag cori-lg-data my.wdl inputs.json cori


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
    jaws run submit my.wdl inputs.json cori


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

When monitoring the runs with :bash:`jaws run status`, each task transitions between the following states. 

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


**Checking the status of a task**

.. code-block:: text

    # the two status commands show the current status of the run or tasks of the run
    jaws run status 121
    jaws run task-status 121

    # the log commands include all the past states of either the run or tasks of the run
    jaws run log 121
    jaws run task-log 121

**Get current or old history of jobs owned by you**

.. code-block:: text

   # get list of your currently running jobs
   jaws run queue                                      
   
   # view history of your jobs for last 7 days 
   jaws run history --days 7


