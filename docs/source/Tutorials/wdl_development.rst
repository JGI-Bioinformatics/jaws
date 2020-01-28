================================
Setting up a Testing Environment
================================

.. role:: bash(code)
   :language: bash

This page describes how to run a workflow from a `cromwell <https://cromwell.readthedocs.io/en/stable/>`_ command, which is the prefered way to test and develop your WDL.
Running your WDL in JAWS will come later.

.. note:: The Workflow Description Language (WDL) 
	WDL is developed by an open source community `openwdl.org <openwdl.org>`_. It is essentially a wrapper around the commands in your pipeline code.  Cromwell is a workflow engine that executes the commands in the WDL.

You can clone a repository that has all the required example code and data.  Then you will run a simple WDL, `align.wdl` that was derived from the `script.sh` file. The steps in the workflow will align some fastq reads to a `reference.fasta` and produce a bam file. 

*******************
Install Miniconda3
*******************
When you need to create an environment for your WDL to run in, docker is a good choice; however, conda environments are good as a preliminary step for testing, and then you can make a docker image easily from the conda env.  

If you don't already have miniconda...

Follow directions at 
https://docs.conda.io/en/latest/miniconda.html


******************************************
Install cromwell or Use Cori Installation
******************************************
Cromwell is a workflow engine that understands and runs WDLs.

You can use the previously installed cromwell at:

.. code-block:: bash

	/global/dna/projectdirs/DSI/workflows/cromwell/java/cromwell.jar


Create a conda environment to install dependencies in

.. code-block:: bash

  conda create -n cromwell python=3
  conda activate cromwell

If you want to install your own cromwell package, then...

.. code-block:: bash

  conda install -y -c bioconda cromwell==0.40

*******************************
Install Other Tool Dependencies
*******************************

.. code-block:: bash

    conda install -y -c bioconda bbmap==38.69
    conda install -y -c bioconda samtools==1.9


***********************************
Download an Example WDL repository
***********************************

.. code-block:: bash

  git clone git@gitlab.com:jfroula/jaws-example-wdl.git
  cd jaws-example-wdl/5min_example

********************
Run the WDL Workflow
********************

.. code-block:: bash
  
  # run with previously installed version
  java -jar /global/dna/projectdirs/DSI/workflows/cromwell/java/cromwell.jar run align.wdl -i inputs.json

  # run with your installed version
  cromwell run align.wdl -i inputs.json


You should see a directory `cromwell-executions`.
The resulting bam file from the alignment is here `cromwell-executions/bbtools/<some-long-hash>/call-samtools/execution/test.sorted.bam`

