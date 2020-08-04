============================================
Setting up a Testing Environment Using Conda
============================================

.. role:: bash(code)
   :language: bash


This page describes how to run a workflow from a `Cromwell <https://Cromwell.readthedocs.io/en/stable/>`_ 
command, which is the prefered way to test and develop your WDL.

Cromwell is a workflow engine that understands and runs WDLs.

You can clone a repository that has all the required example code and data.  Then you will run a simple WDL, `align.wdl` that was derived from the `script.sh` file. The steps in the workflow will align some fastq reads to a `reference.fasta` and produce a bam file. 


.. _install_miniconda:


*******************************************
Install Miniconda and Download Dependencies
*******************************************

* Define all your dependencies and create an environment that your pipeline can be run in. 
* Conda environments are good for testing your code and to ensure it will be portable to a Dockerfile (see `conda overview <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html?highlight=environment>`_).  

How to find packages you need for building environments
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

1. search google "conda install <your_software_name>" as this could give you alternative (private) channels that have your package
2. If that fails, try to search pip packages.

`conda: <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html?highlight=environment>`_ can install python and non-python module (e.g. perl modules or R packages)   

`pip: <https://docs.python.org/3/installing/index.html>`_ installs only python packages but can be used to install stuff in your conda env.  You can search for pip packages at `pypi <https://pypi.org/>`_


Create an environment called bbtools to test our script (script.sh).
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you don't already have conda installed, download `minconda <https://docs.conda.io/en/latest/miniconda.html>`_ and choose the appropriate miniconda.
I chose "64-bit (bash installer)" for MaxOS, Python 2.7.

Installation

.. code-block:: bash

   bash Miniconda2-latest-MacOSX-x86_64.sh
   # follow the directions for installation.  
   
   # Make sure the installation is in your path. 
   # You may need to add <your_installation>/bin to PATH 
   # depending on what option you chose during the installation.
   which conda

*******************************
Install Other Tool Dependencies
*******************************
Create a conda environment to install dependencies in

.. code-block:: bash

  conda create -n test_env python=3
  conda activate test_env  # if this fails, use source activate test_env

Install the dependencies

.. code-block:: bash

    conda install -y -c bioconda bbmap==38.69
    conda install -y -c bioconda samtools==1.9

Of course you can install Cromwell if you wanted to test on your labtop, etc., by adding the command

.. code-block:: bash

    conda install -y -c bioconda Cromwell

***********************************
Download the Example WDL repository
***********************************

.. code-block:: bash

  git clone https://code.jgi.doe.gov/advanced-analysis/jaws-tutorial-examples.git
  cd jaws-tutorial-examples/5min_example

********************
Run the WDL Workflow
********************

.. code-block:: bash
  
  # run with your installed version
  Cromwell run align.wdl -i inputs.json


You should see a directory `Cromwell-executions`.
The resulting bam file from the alignment is here `Cromwell-executions/bbtools/<some-long-hash>/call-samtools/execution/test.sorted.bam`

