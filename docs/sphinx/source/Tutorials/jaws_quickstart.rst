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

Edit `~/jaws.conf` and add the token that you recieved from the JAWS admin.

.. code-block:: text

      [USER]
      token: <your token>

2. Set up the virtual environment.
You will use an existing file. This gives you access to all the JAWS commands.  
By using a symlink, we can update the file without requiring you to re-copy it. 

.. code-block:: text

    ln -s /global/cfs/projectdirs/jaws/jaws-prod/jaws-prod.sh ~

3. Source the environment to activate JAWS commands.

.. code-block:: text

    source ~/jaws-prod.sh

4. To deactivate the environment, use:

.. code-block:: text

    deactivate
    
**************************
Run an Example WDL in JAWS
**************************

Currently JAWS can run on the following resources:

  * CORI (at NERSC)
  * JGI (at LBNL)

When submitting a JAWS run, you must specify the resource to use (i.e. CORI or JGI)

1. Activate the environment you set up above:

.. code-block:: text

    source ~/jaws-prod.sh

2. Clone the example code:

.. code-block:: text

    git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
    cd jaws-tutorial-examples/quickstart

3. List all the sites available to JAWS:

.. code-block:: text

    jaws list-sites  

4. Submit a workflow using JAWS

.. code-block:: text

    jaws submit --tag metagenome_alignment align.wdl inputs.json cori  # note that case doesn't matter for sites. CORI and cori both work.

    # you should see something like this
    100%|███████████████████████████████████| 2929/2929 [00:00<00:00, 1081055.65it/s]
    Copied 2929 bytes in 0.0 seconds.
    100%|███████████████████████████████████| 792/792 [00:00<00:00, 349231.37it/s]
    Copied 792 bytes in 0.0 seconds.
    {
    "max_ram_gb": 5,
    "run_id": 35970
    }

******************
Monitoring the Job
******************

From the output above, we see that the run_id was `35970`.

.. code-block:: text

    # make sure you remember the id of the job submission,
    # if you didn't you can run this to see your run's id
    jaws queue
    
    # check jaws status 35970

    # check status of the tasks (the last command has the most detail)
    jaws task-log 35970
    jaws task-summary 35970
    
***************
Get the results
***************

Once the run status has changed to "download complete", you can write the output to a folder of your choice. This command run without the :bash:`--complete` flag does not copy all the files in the cromwell :bash:`execution` directory, but only the files listed in the :bash:`output{}` section of the WDL. 

.. code-block:: text

    # copy the output of run 1367 to a folder of your choice
    jaws get 35970 my-test-run
    or
    jaws get --complete 35970 my-test-run-complete

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

1. Use the :bash:`errors` command. This should show the contents of the stderr and stdout files created per task from Cromwell. It should only show content when there is an error code >0. 
Sometimes a script will write errors to stdout which will be caught, but sometimes it will correctly write to stderr but return an error code of 0, in which case this command won't show anything.

.. code-block:: text

    jaws errors 35975

2. If there is no error from the errors command or it is not clear, you can manually check the contents of the stderr, stdout, script and script.submit files that are created within each task's working directory (saved in your specified output directory). Following the above example, these files would be in:

.. code-block:: text

    jaws get --complete 35975 my-test-run_error
    my-test-run_error/call-samtools/execution/stderr.submit

The :bash:`script.submit` file is what cromwell used to run the :bash:`script` file.


3. The :bash:`task-log` command can show errors created by the backend (i.e. HTCondor), like timeout errors that occur when your task's runtime section didn't request enough time. 

.. code-block:: text

    jaws task-log 35975
    
    #NAME              CACHED  STATUS  QUEUED               RUNNING              FINISHED             QUEUE_DUR  RUN_DUR  
    bbtools.samtools   False   Failed  2022-06-03 16:22:33  2022-06-03 16:22:34  2022-06-03 16:22:35  0:00:01    0:00:01  
    bbtools.alignment  False   Done    2022-06-03 16:21:23  2022-06-03 16:21:24  2022-06-03 16:22:30  0:00:01    0:01:06

