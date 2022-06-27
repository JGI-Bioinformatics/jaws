============================================
Setting up a Testing Environment Using Conda
============================================

.. role:: bash(code)
   :language: bash

*******
Summary
*******
This section describes how to set up a running environment for developing WDLs using `conda <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html?highlight=environment>`_.  Once the environment is set up, we can run a WDL directly using `Cromwell <https://Cromwell.readthedocs.io/en/stable/>`_ which is a workflow engine that understands WDLs, and is what JAWS uses under the hood. You'll want to develop and test your WDL with Cromwell and only run with JAWS once everything is working.

In this tutorial, you will
	1. build a development environment
	2. run a WDL using cromwell


How to find packages you need for building environments
+++++++++++++++++++++++++++++++++++++++++++++++++++++++
Common ways to install software into docker containers are:

:1. **apt-get**:                 apt-get is useful for installing dependencies and works for ubuntu, which is used for this docker build. For alpine or other operating systems, find the equivalent installer.
:2. **pip**:                     for installing python packages
:3. **conda**:                   installs packages from any language, but the available versions of a particular software may be limited
:4. **installing from source**:  the old fashion way of running make && make install

For this tutorial, conda will be the primary way to install our requirements.  The first step would be to:

1. search for a software package at `anaconda.org <https://anaconda.org/>`_. **bioconda** or **conda-forge** are well trusted channels.
2. If that fails, try to search pip packages. As mentioned, `pip: <https://docs.python.org/3/installing/index.html>`_ installs only python packages but can be used to install stuff in your conda env.  You can search for pip packages at `pypi <https://pypi.org/>`_.
3. If that fails, you may have to install from source.


Install Conda via Miniconda3
++++++++++++++++++++++++++++
You can download `minconda <https://docs.conda.io/en/latest/miniconda.html>`_ and choose the appropriate miniconda.

Installation

.. code-block:: text

   # pick the latest version, for Mac in this case.
   bash Miniconda3-latest-MacOSX-x86_64.sh
   # follow the directions for installation.

   # Choose yes for: "Do you wish the installer to initialize Miniconda3 by running conda init?"
   # This should add the miniconda installation to your PATH by modifying .bash_profile.

   # verify installation
   which conda



Create your Conda Environment
++++++++++++++++++++++++++++++
We will create a conda environment called :bash:`bbtools`

.. code-block:: text

   # create bbtools env
   conda create --name bbtools

   # Now activate your environment.
   conda activate bbtools

Install necessary dependencies into your environment

.. code-block:: text

   # install dependencies
   conda install -c bioconda bbmap==38.84
   conda install -c bioconda samtools==1.11

   # of course you can install Cromwell if you wanted to develop WDLs on your labtop
   conda install -y -c bioconda cromwell

************************
Testing your Environment
************************
Download the Example WDL repository

.. code-block:: text

  git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
  cd jaws-tutorial-examples/5min_example

The following command should run an alignment.

.. code-block:: text

	./script.sh ../data/sample.fastq.bz2 ../data/sample.fasta

This should create file called :bash:`test.sorted.bam`.  We will use this :bash:`script.sh`` to create a WDL later.


Run the WDL Workflow
++++++++++++++++++++
You can run a WDL directly on the command-line (outside of JAWS) by using a Cromwell executable. Either you install your own, i.e using conda, or if you are on CORI, then there is an installation at :bash:`/global/cfs/projectdirs/jaws/cromwell/cromwell.jar` which should point to the latest version used in JAWS.

.. raw:: html

    <details>
    <summary style=color: #448ecf;>Installing your own Cromwell</summary>

.. code-block:: text

    conda activate bbtools # make sure you are in a conda environment first
    conda install cromwell
    cromwell --version

.. raw:: html

    </details>
    <br><br>


**Running with your own conda version**

Make sure the bbtools conda environment is activated and you are in 5min_example.

.. code-block:: text

  # run with your installed version
  cromwell run align.wdl -i inputs.json


**Running with Cori's version**

You can also run the WDL on Cori and you could use the pre-installed version of cromwell.  Make sure the bbtools conda environment is created on cori, is activated and you are in 5min_example.

.. code-block:: text

  # run with your installed version
  java -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run align.wdl -i inputs.json


You should see a directory `Cromwell-executions`.
The resulting bam file from the alignment is here `cromwell-executions/bbtools/<some-long-hash>/call-samtools/execution/test.sorted.bam`
