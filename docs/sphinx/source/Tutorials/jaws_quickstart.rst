===============
JAWS Quickstart
===============

.. role:: bash(code)
  :language: bash

*******
Summary
*******

The following commands assume they are to be run by JGI personell and run on NERSC's Cori cluster.

Set up the environment to start running a pipeline in JAWS.
-----------------------------------------------------------

These are the steps (in order) you'll be required to complete.

.. |check| raw:: html

    <input checked=""  type="checkbox">

.. |br| raw:: html

  <br/>

|check| get a NERSC account |br|
|check| get JAWS token |br|
|check| activate the JAWS virtual environment |br|
|check| get Globus token

***********
Begin Setup
***********

`(see video about setting up the JAWS environment) <https://youtu.be/7qXpMNdQjdw>`_

Before you go further, you need to complete these three steps: 

1) get an account at NERSC.  

    - If you do not have a NERSC account yet, please get one by writing to consult@NERSC.gov

2) get a JAWS token from JAWS admin jaws-support@lbl.gov 

|

*********************************
Activate JAWS Virtual Environment
*********************************


Currently JAWS can run on the following resources:

  * CORI (at NERSC)
  * JGI (at LBNL)

.. note::
    When submitting a JAWS run, you must specify the resource to use (i.e. CORI or JGI)

Do the following

.. code-block:: text

    cp /global/cfs/projectdirs/jaws/jaws-prod/jaws.conf ~/jaws.conf
    chmod 600 ~/jaws.conf

    Edit ~/jaws.conf and add values for the [USER] variables:
      token : This should be the token you got from the JAWS admin

    # Set up the virtual environment
    # You will use an existing one. This gives you access to all the jaws commands.  # By using a symlink, we can update the file without requiring you to re-copy the file.
    ln -s /global/cfs/projectdirs/jaws/jaws-prod/jaws-prod.sh ~

    source ~/jaws-prod.sh
    (use "deactivate" to get out of the environment)

|

***************
Run WDL in JAWS
***************

.. code-block:: text

    # activate the environment you set up above
    source ~/jaws-prod.sh

    # clone the example code
    git clone https://code.jgi.doe.gov/advanced-analysis/jaws-tutorial-examples.git
    (You'll need to use your LBL LDAP credentials).

    cd jaws-tutorial-examples/quickstart

    # run "jaws run submit <workflow> <inputs> <site>"
    jaws run list-sites  # you should see all the sites available to JAWS
    jaws run submit align.wdl inputs.json cori  # note that case doesn't matter for sites.

    # you should see something like this
    2020-11-13 17:51:20,444 - INFO - workflow - Validating WDL, /global/cscratch1/sd/jfroula/JAWS/jaws-tutorial-examples/quickstart/align.wdl
    2020-11-13 17:51:24,762 - INFO - workflow - Maximum RAM requested is 5Gb
    2020-11-13 17:51:24,790 - INFO - workflow - Writing file manifest to .../JAWS-scratch/9cfc798e-2015-4cd8-b1ce-75e56f033ccb.tsv
    2020-11-13 17:51:26,919 - INFO - analysis - Submitted run 1367: {'site_id': 'CORI', 'submission_id': '9cfc798e-2015-4cd8-b1ce-75e56f033ccb', 'input_site_id': 'CORI', 'input_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'output_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'output_dir': '/global/cscratch1/sd/jfroula/JAWS/jaws-tutorial-examples/quickstart/out'}
    {
      "run_id": 1367,
      "site_id": "CORI",
      "status": "uploading",
    } 

******************
Monitoring the Job
******************

From the output above, we see that the run_id was 1367.

.. code-block:: text

    # make sure you remember the id of the job submission,
    # if you didn't you can run this to see your run's id
    jaws run queue
    
    # check jaws status
    jaws run status 1367

    # check status of the tasks (the last command has the most detail)
    jaws run task-status 1367
    jaws run task-log 1367

***************
Get the results
***************

Once the run status has changed to "download complete", you can write the output to a folder of your choice using:

.. code-block:: text

    # copy the output of run 1367 to a folder of your choice
    jaws run get 1367 $SCRATCH/my-test-run


***********
Output
***********

Cromwell will create a directory structure that looks like this: (different from what you'll see):

.. figure:: /Figures/crom-exec.svg
    :scale: 100%

Each task of your workflow gets run inside the :bash:`execution` directory so it is here that you can find any output files including the stderr, stdout & script file. Cromwell is run on scratch and when it is finished, everything below the "cromwell generated hash" is copied to your specified output directory. 

    
So for our theoretical submission

.. code-block:: text

    jaws run submit align.wdl inputs.json cori  

We should see an output folder that looks like this:

.. figure:: /Figures/crom-exec-jaws.svg
    :scale: 100%


Further Debugging Ideas
-----------------------

1) The :bash:`metadata` command will show you the output from the Cromwell server which may have additional debugging information.  Look for "causedBy" message as shown below. This error doesn't tell you much so the next step would be 2) below.

.. code-block:: text

    jaws run metadata 80

    "causedBy": [],
        "message": "Job jgi_dap_leo.assignGenes:4:1 exited with return code 79 which has not been declared as 
        a valid return code. See 'continueOnReturnCode' runtime attribute for more details."
    }

2) Use the :bash:`errors` command. This should show the contents of the stderr file, but only when there was an error code >0. 
Sometimes a script will write to stderr but return an error code of 0, so this command won't show anything.

.. code-block:: text

    jaws run errors 1186


3) Check the contents of the stderr, stdout files that are created within each task's working directory (saved in your specified output directory). Following the above example, your stderr/stdout files would be in:

.. code-block:: text

    out/call-setup/execution/stderr

It is also useful to examine the file called :bash:`script` since this is exactly what cromwell ran.


4) Use the :bash:`task-log` command to show errors that JTM catches, like timeout errors that occur when your task's runtime section didn't request enough time. We are aware of an issue with this command having a long delay, so please be patient until we can re-design the way task-log (and task-status) works.

.. code-block:: text

    jaws run task-log 1186
    
    "jgi_dap_leo.assignGenes 1   5132    running failed  2020-10-28 21:11:14 failed with timeout"

