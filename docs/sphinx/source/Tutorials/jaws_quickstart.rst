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

Currently JAWS can be run at 

  * Cori(NERSC)
  * LBNL(Lawrencium)  

.. note::
	Important  
	When running a JAWS job, you will always start it on Cori, and use the flag --site=[nersc|lbnl] to specify where you want to run it.


So first thing is, you need to have an account at NERSC.

If you do not have a NERSC account yet, please get one by writing to consult@nersc.gov .  

.. code-block:: bash

    # get access to the "activate" command
    module load python  
    
    # TODO: need instructions how to install JAWS
    python3 -m venv ./jawsenv
    source ./jawsenv/bin/activate
    pip install jaws

    # verify you have the command
    jaws

****************************
Set up Globus (do this once)
****************************
You need to set up globus so your files can be transfered between different filesystems (even if you are running everything on cori).  

Setting up globus endpoints

.. code-block:: bash

    1. open globus and create an account if you don't already have one globus.org
    2. once you are logged in, you need to add an endpoint (the institution should have an endpoint already; its where data will be transfered to and from).  
      a. Click on ENDPOINTS in the menu on the left. Note, you should be in the Recently Used tab which is default.  
      b. Search for NERSC DTN   
      c. click on the arrow > at the right of NERSC DTN which is the endpoint details (if you curser over it)   
      d. on the right of the screen, click on Activation (or Extend Activation) You will have to re-activate every 11 days. Follow the directions to authenticate using NERSC credentials.  
    
Setting up globus linked accounts

.. code-block:: bash

    You can link accounts like your NERSC and LBL account. Linking the NERSC account is required for globus to know that its ok to upload and download your data when you are using jaws and thus logged in as NERSC credentials. 
    1. click on ACCOUNT in the left menu.  You should be in the Identities tab. 
    2. click on Link Another Identity
    3. search for NERSC and click continue....follow the authentication steps.  
       a. You should see <yourusername>@nersc.gov.   


************************
Getting JAWS Credentials
************************
The first time you run JAWS, you will have to get a globus token for transfering data (do this once). Run the following command and follow directions.

.. code-block:: bash

    jaws run queue  # this command is to see running jobs but it works to trigger the new user login sequence

You need to follow the directions to get the token, for example:

.. code-block:: bash

    Creating new config file: /global/homes/j/jfroula/.jaws.dev.ini
    Updating user config file at /global/homes/j/jfroula/.jaws.dev.ini
    Welcome, new user
    Do you have a Globus account? [y/n]: y
    Have you been added to the JAWS group? [y/n]: y
    Step 1: Authenticate your NERSC endpoint at https://app.globus.org/file-manager?origin_id=9d6d994a-6d04-11e5-ba46-22000b92c6ec
    Have you activated the nersc#dtn endoint? [y/n]: y
    Step 2: Grant JAWS permissions to transfer files on your behalf
    Please go to this URL and log in
    (HINT: try command-click on the link):
    
    https://auth.globus.org/v2/oauth2/authorize?client_id....

    Then paste the authorization code here:

***************
Run WDL in JAWS
***************

.. code-block:: bash

    # clone the example code
    git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
    
    cd jaws/examples/create_wdl_tutorial
    
    # run jaws run submit <workflow> <inputs> <outdir>
    jaws run submit align.wdl inputs.json out
    
    # you should see something like this
    Submitting to: NERSC
    {
      "run_id": 157
    }
    

******************
Monitoring the Job
******************

.. code-block:: bash

    # make sure you remember the id of the job submission, if you didn't you can run this to see your run's id
    jaws run queue
    
    # check jaws status
    jaws run status 157

***********
Output
***********
All output files should be in "out" in this case.

For debugging, check JAWS run directory.
----------------------------------------

.. code-block:: bash

    # look for "workflowRoot" near bottom of metadata output. 
    # This is the path to where cromwell results reside.
    jaws run metadata 157

