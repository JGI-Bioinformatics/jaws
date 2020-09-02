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
|check| get Globus account and link your NERSC account to Globus |br|
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

3) get an account at globus.org

    - `for help logging in <https://docs.globus.org/how-to/get-started>`_

*************
Set up Globus 
*************

You need to set up Globus so your files can be transfered between different filesystems. 

Setting up Globus endpoints

.. code-block:: bash

    1. open globus.org (you should have created an account already for the previous steps to work).

      a. under the "Account" tab, make sure you have something like <username>@globusid.org

    2. once you are logged in, you need to add an endpoint (institutions that JAWS
       uses have endpoints already; its where data will be transfered to and from).

      a. Click on "ENDPOINTS" in the menu on the left. Note, you should be in the
         "Recently Used" tab which is default.

      b. Search for NERSC DTN

      c. Click on the arrow > at the right of NERSC DTN which is the endpoint details 
         (if you curser over it)

      d. On the right of the screen, click on "Activation" (or "Extend Activation")
         Follow the directions to authenticate using NERSC credentials!!
         You will have to re-activate every 11 days.


Setting up Globus linked accounts 

.. code-block:: bash

    You can link accounts like your NERSC and LBL account. 
    Linking the NERSC account is required for Globus to know that its ok to upload 
    and download your data when you are using jaws and thus logged in as NERSC credentials. 

    1. Click on "ACCOUNT" in the left menu.  You should be in the "Identities" tab. 

    2. Click on "Link Another Identity"

    3. Search for NERSC and click continue....follow the authentication steps.  

       a. You should see <yourusername>@NERSC.gov listed.   

.. warning:: 
	You need to re-activate your Globus Endpoint every 11 days.  JAWS should give you an appropriate error if you need to do this. Go to globus.org and click on "ENDPOINTS".  If "NERSC DTN" says "inactive", you can click on the activate endpoint symbol at the right.

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

.. code-block:: bash

    cp /global/cfs/projectdirs/jaws/jaws-prod/jaws.conf ~
    chmod 600 ~/jaws.conf

    Edit ~/jaws.conf and add values for the [USER] variables:
      token : This should be the token you got from the JAWS admin
      staging_dir : Set to a JAWS subdir in your scratch dir, e.g. /global/cscratch1/sd/YOURUID/jaws

    # Set up the virtual environment
    # You will use an existing one. This gives you access to all the jaws commands.
    # By using a symlink, we can update the file without forcing you to re-copy the file.
    ln -s /global/cfs/projectdirs/jaws/jaws-prod/jaws-prod.sh ~

    source ~/jaws-prod.sh
    (use "deactivate" to get out of the environment)

    # Get the jaws-auth token. 
    # After running this command, follow directions to get a token from Globus.
    jaws login

|

***************
Run WDL in JAWS
***************

.. code-block:: bash

    # activate the environment you set up above
    source ~/jaws

    # clone the example code
    git clone https://code.jgi.doe.gov/advanced-analysis/jaws-tutorial-examples.git

    cd jaws-tutorial-examples/quickstart

    # run "jaws run submit <workflow> <inputs> <full path to outdir> <site>"
    jaws run list-sites  # you should see CORI
    jaws run submit align.wdl inputs.json out cori  # note that case doesn't matter for sites.

    # you should see something like this
    2020-04-16 13:04:18,434 - INFO - workflow - Validating WDL, align.wdl
    2020-04-16 13:04:20,357 - INFO - workflow - Validating inputs file, inputs.json
    2020-04-16 13:04:22,084 - INFO - workflow - Maximum RAM requested is 0Gb
    2020-04-16 13:04:22,085 - INFO - workflow - Staging WDLs to <fullpath>/JAWS-scratch
    2020-04-16 13:04:22,088 - INFO - workflow - Staging infiles to <fullpath>/JAWS-scratch/CORI
    2020-04-16 13:04:22,093 - INFO - workflow - Writing file manifest to <fullpath>/JAWS-scratch/ca626c3e-ad65-44b8-a55a-4ce310d2108b.tsv

    {
        "output_dir": "<fullpath>/examples/create_wdl_tutorial/out",
        "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
        "run_id": 80,
        "site_id": "CORI",
        "status": "uploading",
        "submission_id": "ca626c3e-ad65-44b8-a55a-4ce310d2108b",
        "upload_task_id": "77810d8e-801d-11ea-97a5-0e56c063f437"
    }
    

******************
Monitoring the Job
******************

From the output above, we see that the run_id was 80.

.. code-block:: bash

    # make sure you remember the id of the job submission,
    # if you didn't you can run this to see your run's id
    jaws run queue
    
    # check jaws status
    jaws run status 80

    # check status of the tasks (the last command has the most detail)
    jaws run task-status 80
    jaws run task-log 80


***********
Output
***********

All output files and logs should be in the output directory that you specified, "out" in this case.


If a Job Fails
--------------

If a job fails, your output dir will contain a copy of the raw `Cromwell <https://cromwell.readthedocs.io/en/stable/>`_ execution directory. 

For example, a directory like this will exist:

.. figure:: /Figures/crom-exec.svg
    :scale: 100%

You will have to look at the task's stdout, stderr & script files to see what went wrong.

Further Debugging Ideas
-----------------------

.. code-block:: bash

    # The output command should show you the contents of the stderr, stdout (same content as the stderr mentioned above).
    jaws run output 80

    # The metadata command will show you the output from the Cromwell server which may have additional debugging information.
    # look for "causedBy"
    jaws run metadata 80
