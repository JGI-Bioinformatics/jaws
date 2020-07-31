======================== 
How to Run a WDL in JAWS
========================

.. role:: bash(code)
   :language: bash

.. note::
   Before these commands will work, you need to set up everything:
   See `Quickstart Example </Tutorials/jaws_quickstart.html>`_ for full setup instructions.


Now you can submit a run

:bash:`jaws run submit <wdl> <inputs json> <full path to out dir> <nersc|lbnl>` 


*************
JAWS commands
*************


There are four top level JAWS commands (but many sub commands):

  :bash:`jaws login` after running this command, you need to follow the directions to obtain a globus token.

  :bash:`jaws run` command for job management (e.g. submission, monitoring and log-viewing). 

  :bash:`jaws wdl`  command for workflow management (e.g. add, delete & view existing workflows). 

  :bash:`jaws status`  displays the status of some services at different sites (i.e. nersc & lbnl).


jaws *login*:
----------------------

You only need to run this once in theory to get a globus token.  Follow the directions to obtain the token from globus and this will allow JAWS to transfer data around via globus in your name.

.. code-block:: bash

    jaws login
    
    # This will generate a URL that you need to follow to get the globus token 
    


jaws *run* options:
-------------------

.. code-block:: bash

   cancel <run_id>                cancel a run, prints whether aborting was successful or not
   output <run_id>                view the stdout/stderr for all tasks
   errors <run_id>                view the stdout/stderr for failed tasks
   history <num_days>             list past runs for a given number of days including results
   list-sites                     list available sites you can run jobs from
   log                            view the history of Run state transitions
   metadata <run_id>              detailed information about a run
   queue                          list your unfinished runs
   status <run_id>                show current status of a single run including result (succeeded/failed)
   submit <wdl> <inputs>          submit a run for execution
   task-status <run_id>           show current status of each task of a run
   task-log                       view the history of tasks' state transitions

   
jaws *wdl* options:
-------------------

.. code-block:: bash

  about       return README document for a workflow
  add         add a workflow to the catalog. do this before **release**
  delete      remove a workflow from the catalog
  get         get WDL specification for a workflow
  list        list shared workflows
  release     mark a version as immutable production release
  update-doc  update a workflow README in the catalog
  update-wdl  update a workflow WDL in the catalog
  versions    list available versions of a workflow
 
 
Note: The WDL of 'released' workflows cannot be changed.


jaws *status*:
----------------------

This command shows the status of the different services. Some services are site specific; if a services at LBNL is down, running JAWS at the NERSC site should still work.

.. code-block:: bash

    {
      "JAWS-Central": "UP",
      "LBNL-Cromwell": "Unknown",
      "LBNL-RMQ": "UP",
      "LBNL-Site": "DOWN",
      "NERSC-Cromwell": "UP",
      "NERSC-RMQ": "UP",
      "NERSC-Site": "UP"
    }



Examples
--------

**Anyone can share a WDL. To see a list of workflows available in the catalog run**

::

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

::

   # note that a version is required
   jaws wdl about fq_count dev 


**To run a wdl**

::

    # submit it
    jaws run submit my.wdl inputs.json out nersc


**To run a WDL from the catalog, there are a couple extra steps (from "jaws run list" we saw there is a wdl in the catalog called fq_count)**

::

    # create the wdl
    jaws wdl get fq_count dev > my.wdl
    
    # create a template for inputs.json 
    jaws run inputs fq_count > inputs.json

    # cusomize the values in inputs.json
    vi inputs.json

    # run as usual
    jaws run submit my.wdl inputs.json out nersc


.. note::

    From any job submition, you can see a run id (i.e. below you can see 121). Use this for future commands.

::

  # output looks like
  {
  "output_dir": "<full_path>/out",
  "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
  "run_id": 121,
  "site_id": "NERSC",
  "status": "uploading",
  "submission_id": "7d2606b9-569f-4d50-9423-c1acb5441c6b",
  "upload_task_id": "07ffa460-88ac-11ea-b3ba-0ae144191ee3"
  }

 


**See the status & metadata of a run using job ID**

::

  jaws run status 121

  # there's some usefull stuff in metadata so check it out
  jaws run metadata 121


**Get current or old history of jobs owned by you**

::

   # get list of your currently running jobs
   jaws run queue                                      
   
   # view history of your jobs for last 7 days 
   jaws run history --days 7
