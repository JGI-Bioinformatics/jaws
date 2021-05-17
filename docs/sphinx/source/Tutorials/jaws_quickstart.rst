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

|check| get JAWS token |br|
|check| symlink the file responsible for starting the virtual environment |br|
|check| copy the configuration file that will contain your token. |br|
|check| activate the JAWS virtual environment |br|

***********
Begin Setup
***********

.. note::
    Before you can complete the JAWS setup, you need to get a **JAWS token** from JAWS admin jaws-support@lbl.gov 

***********
Set up JAWS
***********

1. Set up the config file that  has the JAWS token in it.

.. code-block:: text

    cp /global/cfs/projectdirs/jaws/jaws-prod/jaws.conf ~/jaws.conf
    chmod 600 ~/jaws.conf

Edit ~/jaws.conf and add the token that you recieved from the JAWS admin.

.. code-block:: text

      [USER]
      token: <your token>

2. Set up the virtual environment.
You will use an existing file. This gives you access to all the jaws commands.  
By using a symlink, we can update the file without requiring you to re-copy it. 

.. code-block:: text

    ln -s /global/cfs/projectdirs/jaws/jaws-prod/jaws-prod.sh ~

3. Source the environment to activate JAWS commands

.. code-block:: text

    source ~/jaws-prod.sh
    (use "deactivate" to get out of the environment)

|

**************************
Run an Example WDL in JAWS
**************************

Currently JAWS can run on the following resources:

  * CORI (at NERSC)
  * JGI (at LBNL)

When submitting a JAWS run, you must specify the resource to use (i.e. CORI or JGI)

.. code-block:: text

    # activate the environment you set up above
    source ~/jaws-prod.sh

    # clone the example code
    git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git

    cd jaws-tutorial-examples/quickstart

    # you should see all the sites available to JAWS
    jaws list-sites  

    jaws submit --tag metagenome_alignment align.wdl inputs.json cori  # note that case doesn't matter for sites. CORI and cori both work.

    # you should see something like this
    2020-11-13 17:51:20,444 - INFO - workflow - Validating WDL, /global/cscratch1/sd/jfroula/JAWS/jaws-tutorial-examples/quickstart/align.wdl
    2020-11-13 17:51:24,762 - INFO - workflow - Maximum RAM requested is 5Gb
    2020-11-13 17:51:24,790 - INFO - workflow - Writing file manifest to .../JAWS-scratch/9cfc798e-2015-4cd8-b1ce-75e56f033ccb.tsv
    2020-11-13 17:51:26,919 - INFO - analysis - Submitted run 1367: {'site_id': 'CORI', 'submission_id': '9cfc798e-2015-4cd8-b1ce-75e56f033ccb', 'input_site_id': 'CORI', 'input_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'output_endpoint': '9d6d994a-6d04-11e5-ba46-22000b92c6ec', 'output_dir': '/global/cscratch1/sd/jfroula/JAWS/jaws-tutorial-examples/quickstart/out'}
    {
      "run_id": 1367,
      "site_id": "CORI",
      "status": "uploading",
      "tag": "metagenome_alignment"
    } 

******************
Monitoring the Job
******************

From the output above, we see that the run_id was 1367.

.. code-block:: text

    # make sure you remember the id of the job submission,
    # if you didn't you can run this to see your run's id
    jaws queue
    
    # check jaws status
    jaws status 1367

    # check status of the tasks (the last command has the most detail)
    jaws task-status 1367
    jaws task-log 1367

***************
Get the results
***************

Once the run status has changed to "download complete", you can write the output to a folder of your choice. This command does not copy the "input" files.

.. code-block:: text

    # copy the output of run 1367 to a folder of your choice
    jaws get 1367 $SCRATCH/my-test-run

Alternatively, this command will display the directory where JAWS has saved the results. It represents all the files from Cromwell, including the "input" files.

.. code-block:: text

    jaws status --verbose 1367

    or

    jaws history


******************************
The Output Directory Explained
******************************

Cromwell will create a directory structure similar to this...

.. figure:: /Figures/crom-exec.svg
    :scale: 100%

Each task of your workflow gets run inside the :bash:`execution` directory so it is here that you can find any output files including the stderr, stdout & script file. Cromwell is run on scratch and when it is finished, everything below the "cromwell generated hash" is copied to your specified output directory. 

    
So for our theoretical submission

.. code-block:: text

    jaws submit align.wdl inputs.json cori  

We should see an output folder that looks like this:

.. figure:: /Figures/crom-exec-jaws.svg
    :scale: 100%


Further Debugging Ideas
-----------------------

1) Use the :bash:`errors` command. This should show the contents of the stderr and stdout files created per task from Cromwell. It should only show content when there is an error code >0. 
Sometimes a script will write errors to stdout which will be caught, but sometimes it will correctly write to stderr but return an error code of 0, in which case this command won't show anything.

.. code-block:: text

    jaws errors 1186


2) If there is no error from the errors command or it is not clear, you can manually check the contents of the stderr, stdout, script and script.submit files that are created within each task's working directory (saved in your specified output directory). Following the above example, these files would be in:

.. code-block:: text

    out/call-setup/execution/stderr

The :bash:`script.submit` file is what cromwell used to run the :bash:`script` file.


3) The :bash:`task-log` command can show errors created by the backend (i.e. JTM), like timeout errors that occur when your task's runtime section didn't request enough time. 

.. code-block:: text

    jaws task-log 1186
    
    "jgi_dap_leo.assignGenes 1   5132    running failed  2020-10-28 21:11:14 failed with timeout"

