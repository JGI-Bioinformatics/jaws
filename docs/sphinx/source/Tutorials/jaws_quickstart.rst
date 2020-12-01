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

.. code-block:: text

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
         You will have to re-activate every 11 days, as is NERSC's policy.

      e. Presently, if you want to use the "JGI" site, you need to add the endpoint "lbnl#lrc".  
         Start from step 2.b. and add "lbnl#lrc".  You'll have to have a Lawrencium account from 
         LAB IT, which is different from your LBNL credentials. This will not be a requirement 
         in future versions. Also, for now, you'll have to re-activate every 7 days, as is LBNL's policy.

After you have added nersc and lbl endpoints, you may not see them listed on your "Endpoints" page under the tab "Recently Used" until you submit a JAWS run.  You can type lbnl#lrc or nersc dtn in the uppermost search window (the one for Endpoints) and they should show up with a "STATUS" of ready. 


Setting up Globus linked accounts 

.. code-block:: text

    You can link accounts like your NERSC and LBL account. 
    Linking the NERSC account is required for Globus to know that its ok to upload 
    and download your data when you are using jaws and thus logged in as NERSC credentials. 

    1. Click on "ACCOUNT" in the left menu.  You should be in the "Identities" tab. 

    2. Click on "Link Another Identity"

    3. Search for NERSC and click continue....follow the authentication steps.  

       a. You should see <yourusername>@NERSC.gov listed.   

.. warning:: 
    You need to re-activate your nersc#dtn Globus Endpoint every 11 days and lbnl#lrc every 7 days. We are working towards removing these requirements.  JAWS should give you an appropriate error if you need to re-activate your token. Go to globus.org and click on "ENDPOINTS".  If "NERSC DTN" says "inactive", you can click on the activate endpoint symbol at the right.

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
      staging_dir : Set to a JAWS subdir in your scratch dir, e.g. /global/cscratch1/sd/YOURUID/jaws

    # Set up the virtual environment
    # You will use an existing one. This gives you access to all the jaws commands.  # By using a symlink, we can update the file without requiring you to re-copy the file.
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

.. code-block:: text

    # activate the environment you set up above
    source ~/jaws-prod.sh

    # clone the example code
    git clone https://code.jgi.doe.gov/advanced-analysis/jaws-tutorial-examples.git
    (You'll need to use your LBL LDAP credentials).

    cd jaws-tutorial-examples/quickstart

    # run "jaws run submit <workflow> <inputs> <full path to outdir> <site>"
    jaws run list-sites  # you should see all the sites available to JAWS
    jaws run submit align.wdl inputs.json out cori  # note that case doesn't matter for sites.

    # you should see something like this
    2020-11-13 17:51:20,444 - INFO - workflow - Validating WDL, /global/cscratch1/sd/jfroula/JAWS/jaws-tutorial-examples/quickstart/align.wdl
    2020-11-13 17:51:24,762 - INFO - workflow - Maximum RAM requested is 5Gb
    2020-11-13 17:51:24,790 - INFO - workflow - Writing file manifest to .../JAWS-scratch/9cfc798e-2015-4cd8-b1ce-75e56f033ccb.tsv
    2020-11-13 17:51:26,919 - INFO - analysis - Submitted run 1367: {'site_id': 'CORI', 'submission_id': '9cfc798e-2015-4cd8-b1ce-75e56f033ccb', 'input_site_id': 'CORI', 'input_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'output_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'output_dir': '/global/cscratch1/sd/jfroula/JAWS/jaws-tutorial-examples/quickstart/out'}
    {
      "output_dir": ".../jaws-tutorial-examples/quickstart/out",
      "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
      "run_id": 1367,
      "site_id": "CORI",
      "status": "uploading",
      "submission_id": "9cfc798e-2015-4cd8-b1ce-75e56f033ccb",
      "upload_task_id": "e8048078-261b-11eb-8fbe-0a34088e79f9"
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


***********
Output
***********

Cromwell will create a directory structure that looks like this: (different from what you'll see):

.. figure:: /Figures/crom-exec.svg
    :scale: 100%

Each task of your workflow gets run inside the :bash:`execution` directory so it is here that you can find any output files including the stderr, stdout & script file. Cromwell is run on scratch and when it is finished, everything below the "cromwell generated hash" is copied to your specified output directory. 

    
So for our theoretical submission

.. code-block:: text

    jaws run submit align.wdl inputs.json out cori  

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

