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

  * Cori and 
  * LBNL (Lawrencium). 

You need to have an account at one of the facilities.

  If you do not have a NERSC account yet, please get one by writing to consult@nersc.gov .  

  If you do not have a LBNL account yet, please get one by following directions at `http://scs.lbl.gov/getting-an-account <http://scs.lbl.gov/getting-an-account>`_

This example is for running on Cori but the steps are the same; except at lbl, you'll need to install miniconda so you have access to the :bash:`activate` command. You can see :ref:`install_miniconda`.

.. code-block:: bash

    # get access to the "activate" command
    module load python  
    
    # create this directory if you don't already have it
    mkdir ~/.conda/envs
    
    # create a local conda environment for jaws to make life easier in the future
    ln -s /global/cfs/projectdirs/jaws/prod/cli/ ~/.conda/envs/jaws
    source activate jaws 
    ... or use this instead if you have conda set up this way (you'll know if you do).
    conda activate jaws

    # verify you have the command
    jaws

*********************
Login (do this once)
*********************

.. code-block:: bash

	jaws login <NERSC password + MFA>

***************
Run WDL in JAWS
***************

.. code-block:: bash

    # clone the example code
    git clone https://gitlab.com/jfroula/jaws-example-wdl.git
    
    cd jaws-example-wdl/create_wdl_tutorial
    
    # run jaws submit <workflow> <inputs> <subworkflow>
    jaws submit align.wdl inputs.json
    
    # you should see something like this
    Successfully queued job c6d403a0-bc07-495f-8250-7132593e7d7d
    

******************
Monitoring the Job
******************

.. code-block:: bash

    # make sure you copy the id of the job submission, if you didn't you can run this to see your run's id
    jaws queue
    
    # check jaws status
    jaws status c6d403a0-bc07-495f-8250-7132593e7d7d

***********
Find Output
***********

Check JAWS results directory.
------------------------------

.. code-block:: bash

    # look for "workflowRoot" near bottom of metadata output. 
    # This is the path to where cromwell results reside.
    jaws metadata c6d403a0-bc07-495f-8250-7132593e7d7d

    # or this will just give you the output files
    jaws output c6d403a0-bc07-495f-8250-7132593e7d7d
    
