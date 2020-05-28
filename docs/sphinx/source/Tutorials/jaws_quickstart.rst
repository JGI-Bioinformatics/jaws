===============
JAWS Quickstart
===============

.. role:: bash(code)
  :language: bash

*******
Summary
*******

This is an example of how you can set up the environment to start running a pipeline in JAWS.

.. note:: 
    Some Definitions:

    * The Workflow Description Language (WDL) is developed by an open source community `openwdl.org <openwdl.org>`_. It is essentially a wrapper around the commands in your pipeline code.  
    * Cromwell is a workflow engine that executes the commands in the WDL.

*******************************
Set up JAWS Environment 
*******************************


`(see video about settting up the JAWS environment) <https://youtu.be/7qXpMNdQjdw>`_

Currently JAWS can be run at 

  * Cori(NERSC)
  * LBNL(Lawrencium)  

.. note::
    When running a JAWS job, you will always start it on Cori. You can specify as a command argument which site you want [nersc|lbnl].


To begin, you should do these three things 

1) get an account at NERSC.  

	- If you do not have a NERSC account yet, please get one by writing to consult@nersc.gov .  

2) get an account at globus.org

	- `for help logging in <https://docs.globus.org/how-to/get-started>`_

3) get a JAWS token from one of the JAWS admins (jlfroula@lbl.gov or eskirton@lbl.gov)


On Cori, do the following

.. code-block:: bash

    rm ~/.jaws-dev.ini  # delete old if exists
    cp /global/cfs/projectdirs/jaws/jaws-dev/jaws.conf ~
    chmod 600 ~/jaws.conf

    Edit ~/jaws.conf and add values for the [USER] variables:
      token : This should be the token you got from the JAWS admin
      staging_dir : Set to a JAWS subdir in your scratch dir, e.g. /global/cscratch1/sd/YOURUID/jaws

    # Set up the virtual environment
    # You will use an existing one. This gives you access to all the jaws commands.
    ln -s /global/cfs/projectdirs/jaws/jaws-dev/ ~

    source ~/jaws-dev/bin/activate
    (use "deactivate" to get out of the environment)

    # Get the jaws-auth token. 
    # After running this command, follow directions to get a token from globus.
    jaws login


*************
Set up Globus 
*************

You need to set up globus so your files can be transfered between different filesystems (even if you are running everything on cori).  

Setting up globus endpoints

.. code-block:: bash

    1. open globus.org (you should have created an account already for the previous steps to work).

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


Setting up globus linked accounts 

.. code-block:: bash

    You can link accounts like your NERSC and LBL account. 
    Linking the NERSC account is required for globus to know that its ok to upload 
    and download your data when you are using jaws and thus logged in as NERSC credentials. 

    1. Click on "ACCOUNT" in the left menu.  You should be in the "Identities" tab. 

    2. Click on "Link Another Identity"

    3. Search for NERSC and click continue....follow the authentication steps.  

       a. You should see <yourusername>@nersc.gov listed.   


***************
Run WDL in JAWS
***************

.. code-block:: bash

    # clone the example code
    git clone https://gitlab.com/jfroula/jaws-quickstart-example.git

    cd jaws-quickstart-example

    # run jaws run submit <workflow> <inputs> <full path to outdir> <site: [nersc|lbnl]>
    jaws run submit align.wdl inputs.json out nersc

    # you should see something like this
    2020-04-16 13:04:18,434 - INFO - workflow - Validating WDL, align.wdl
    2020-04-16 13:04:20,357 - INFO - workflow - Validating inputs file, inputs.json
    2020-04-16 13:04:22,084 - INFO - workflow - Maximum RAM requested is 0Gb
    2020-04-16 13:04:22,085 - INFO - workflow - Staging WDLs to <fullpath>/JAWS-scratch
    2020-04-16 13:04:22,088 - INFO - workflow - Staging infiles to <fullpath>/JAWS-scratch/NERSC
    2020-04-16 13:04:22,093 - INFO - workflow - Writing file manifest to <fullpath>/JAWS-scratch/ca626c3e-ad65-44b8-a55a-4ce310d2108b.tsv

    {
        "output_dir": "<fullpath>/examples/create_wdl_tutorial/out",
        "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
        "run_id": 80,
        "site_id": "NERSC",
        "status": "uploading",
        "submission_id": "ca626c3e-ad65-44b8-a55a-4ce310d2108b",
        "upload_task_id": "77810d8e-801d-11ea-97a5-0e56c063f437"
    }
    

******************
Monitoring the Job
******************

From the output above, we see that the run_id was 80.

.. code-block:: bash

    # make sure you remember the id of the job submission, if you didn't you can run this to see your run's id
    jaws run queue
    
    # check jaws status
    jaws run status 80

***********
Output
***********
All output files and logs should be in "out" in this case.

For debugging
-------------

.. code-block:: bash

    # This command should show you the contents of the stderr, stdout, and 
    # script files created by your task commands
    jaws run errors 80

