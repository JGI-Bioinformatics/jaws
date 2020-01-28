==================================
How to Create and Share Workflows 
==================================

.. role:: bash(code)
   :language: bash


Below is a step-by-step tutorial showing you how to get a workflow from a bash script into a WDL and then into JAWS.
Here is a link that just talks about :doc:`building wdls </Tutorials/building_wdls>`

.. note:: 
   If you intend on adding your workflow(WDL) to JAWS, you will need permissions to https://bitbucket.com/berkeleylab/jgi-workflows so you can clone the repository and push your project. Please contact jlfroula@lbl.gov.


*******
Summary
*******

  * JAWS uses a workflow engine (Cromwell), which requires a pipeline be wrapped in a workflow language
  * We use Workflow Definition Language 

If you take the time to learn WDL, you will be able to add your own piplines to JAWS as well as adapt other published pipelines for your use. Any pipeline registered in JAWS is available outside the JGI.  Anyone will be able to download and run the open source JAWS on their own computers and are free to make customizations.


*************************************************
Things to Consider When Converting your Workflow
*************************************************

Cromwell runs a series of "tasks" that can be handled differently from each other. For instance, some tasks can be run in parallel, run on large memory machines, run under different python versions, etc.

You will want to organize each main step of your pipeline into a self contained wrapper that can run and be tested independently. Each step, or wrapper script will represent a "task" in the WDL. This is not the only way to create a WDL file but will make linking together tasks, and debugging, easier.

Here are the steps we're going to take for this tutorial:

   1. start from a pipeline written in bash 
   2. make the docker image 
   3. create the actual WDL file
   4. register it in JAWS.

For criteria on how you should organize your pipeline code, see the bottom of this page. 

********************************
How to create a workflow in JAWS
********************************
For this tutorial, I will be using the example code from `jaws-example-wdl <https://gitlab.com/jfroula/jaws-example-wdl>`_.
To follow along, do...

.. code-block:: bash

   git clone git@gitlab.com:jfroula/jaws-example-wdl.git
   cd jaws-example-wdl/create_wdl_tutorial
   


1) Create your environment
--------------------------

* Define all your dependencies and create an environment that your pipeline can be run in. 
* You will want to create an environment that can be converted into a Docker image ( `docker tutorial <https://docs.docker.com/get-started/>`_. )
* Conda environments are good for testing your code and to ensure it will be portable to a Dockerfile (see `conda overview <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html?highlight=environment>`_).  

Find packages you need for building environments by 
+++++++++++++++++++++++++++++++++++++++++++++++++++

1. search google "conda <your_software_name>" as this could give you alternative (private) channels that have your package
2. search pip packages.

`conda: <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html?highlight=environment>`_ can install python and non-python module (e.g. perl modules or R packages)   

`pip: <https://docs.python.org/3/installing/index.html>`_ installs only python packages but can be used to install stuff in your conda env.  You can search for pip packages at `pypi <https://pypi.org/>`_

Lets create a conda environment called bbtools to test our script (script.sh).
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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


Create the conda environment

.. code-block:: bash

   # create bbtools env
   conda create --name bbtools

   # Now activate your environment.
   source activate bbtools

Install stuff into your environment

.. code-block:: bash

   # install dependencies
   conda install -y -c bioconda bbmap==38.49
   conda install -y -c bioconda samtools==1.9

3) Testing your environment
---------------------------

For our tutorial, we only have a sample wrapper that will become a "task" in the WDL, script.sh. 

.. note :: 
   Each wrapper should write its output to the **current working directory**. You can copy files to other directories after the pipeline has finished.


Try running the script to test environment. You need to be in this directory <your_repo_clone>/jaws-example-wdl/create_wdl_tutorial/. Also, make sure you activated the bbtools environment.

.. code-block:: bash
   
   ./script.sh reads.fq reference.fasta

This should create a bam file (test.sorted.bam).

4) Create docker image
----------------------

   Next we'll describe how to create a Dockerfile and register it with hub.docker.com. (you'll have to create a repository on `hub.docker.com <hub.docker.com>`_ first).  Follow this link if you need more information on how to `building dockerfiles <https://docs.docker.com/get-started/part2/#define-a-container-with-dockerfile>`_.

   To make the Dockerfile, you can use the same commands you used for the conda environment.  Notice that it is good practice to specify the versions for each software like I have done in the Dockerfile. There may be different versions of a conda package for different operating systems, so don't assume the versions I used will work for your operating system. Of course, you can drop the versions altogether to get the latest version.

The Dockerfile looks like

.. code-block:: bash

   FROM continuumio/miniconda2

   # install software
   RUN conda install -c bioconda bbmap
   RUN conda install -c bioconda samtools

   # this will give us a workingdir within the container (e.g. a place we can mount data to)
   WORKDIR /bbmap

   # move script into container
   COPY script.sh /usr/local/bin/script.sh

Build the image and upload to hub.docker.com. You need to use your docker hub user name to tag the image when you are building it.

.. code-block:: bash

   # create a "Build" directory and create docker container from there so its a small image. Its good practice to always create an image in 
   # a directory containing only the required files.
   mkdir Build 
   cp script.sh Dockerfile Build/
   cd Build
   docker build --tag <your_docker_hub_user_name>/bbtools:1.0.0 .
   cd ../


Test that the script runs in the docker container

.. code-block:: bash

   docker run jfroula/bbtools:1.0.0 script.sh
 
   # if you are in the directory where the data is, this should produce a bam file
   docker run --volume="$(pwd):/bbmap" jfroula/bbtools:1.0.0 script.sh reads.fq reference.fasta


When you are convinced the docker image is good, you can register it with `hub.docker.com <hub.docker.com>`_  (you need to make an account first).

.. code-block:: bash

   docker login
   docker push <your_docker_hub_user_name>/bbtools:1.0.0


5) Test your image on cori
--------------------------

Test the docker container on cori.nersc.gov. You'll need to use shifter instead of docker to run your workflow.

example:

.. code-block:: bash

   # pull image from hub.docker.com
   shifterimg pull jfroula/bbtools:1.0.0

   # run your wrapper script. notice we are running the script.sh that was saved inside the image
   shifter --image=jfroula/bbtools:1.0.0 script.sh


6) Compose the actual WDL
-------------------------

This subject is a tutorial in itself. I will continue with our :bash:`script.sh` example at this link :doc:`building wdls </Tutorials/building_wdls>`

When creating real WDLs, this step should be relatively easy if you were able to neatly isolate each step in the workflow as a wrapper.

Learning WDL syntax is best done through the official `WDL docs <https://software.broadinstitute.org/wdl/documentation/>`_ . Here are some `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.


The script.sh that is supplied with the repo has two essential commands: 

.. code-block:: bash
 
   	# align reads to reference contigs
	bbmap.sh in=$READS ref=$REF out=test.sam

	# create a bam file from alignment
	samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam

And it has two inputs :bash:`READS` and :bash:`REF`

So our WDL should look like this. Remember to replace my docker image name with yours.

.. code-block:: bash

   workflow bbtools {
     File reads
     File ref

     call alignment {
       input: fastq=reads,
              fasta=ref
     }
     call samtools {
       input: sam=alignment.sam
    }
   }

   task alignment {
     File fastq
     File fasta

     command {
        shifterimg pull jfroula/bbtools:1.2.1 && \
        shifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=${fastq} ref=${fasta} out=test.sam
     }
     output {
       File sam = "test.sam"
     }
   }


   task samtools {
     File sam

     command {
       shifter --image=jfroula/bbtools:1.2.1.samtools view -b -F0x4 ${sam} | shifter --image=jfroula/bbtools:1.2.1.samtools sort - > test.sorted.bam
     }
     output {
       File bam = "test.sorted.bam"
     }
   }


For a description of what each section of the code does, please refer to the WDL references I listed above.

This sample WDL is also in the repository, called v1.0.1.wdl



7) Adding your WDL to JAWS
--------------------------

You should already have permissions to clone and push to https://bitbucket.com/berkeleylab/jgi-workflows .
(Make your requests by sending an email to jlfroula@lbl.gov if you don't).

.. code-block:: bash

   git clone git@bitbucket.org:berkeleylab/jgi-workflows.git
   cd jgi-workflows
   
   # Create a folder corresponding to the name of your workflow. This will be the public workflow name.
   mkdir <my_public_workflow_name>

Add at least the following two files:

   (a) a README.md file describing the workflow. This will be the public face of your workflow. Release notes should be added to this file.
   (b) the WDL file named by it's version.  (e.g. :bash:`v2.1.9.wdl`).  You may have multiple WDL files, each corresponding to a different version.  


If you have two versions of your WDL, you might see something like this in your directory

.. code-block:: bash
   
   2.1.0.wdl  2.1.9.wdl	README.md


.. note::
   When you create your own README.md, keep in mind this will be public and provide users with info on what to expect from the workflow. For ideas on some things to include, see this template :doc:`README.md </Tutorials/suggested_readme_template>`

   You can see the README for each workflow by the command :bash:`wf about bbtools/1.0.0`


Now push your files to the main repo

.. code-block:: bash

   git add .
   git commit -m 'added my_public_workflow_name'
   git push


Use this command to see if your wdl was registered with JAWS (may take up to 1 hr. to be added to registry)

.. code-block:: bash

   wf list 


*******************
helpful references
*******************

* :doc:`Hello world type example WDLs </Tutorials/example_wdls>`
* `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.
* :doc:`Best Practices when Writing WDLs </Intro/best_practices>`

How to organize your pipeline code.
-------------------------------

   For each major step in your workflow, create a wrapper script that will become a "task" in the WDL. 
   
   In deciding how to split your workflow, consider:

     a. group code that will have the same compute requirements(e.g. large memory machine).
     b. group code that can be run in parallel.
     c. does a group of code make it possible to have a simple input and output structure (e.g. one file in, one file out).
     d. does a group of code make a sub-wdl that you can re-use in other workflows (e.g. alignment or assembly).
     e. does a group of code make it easier to test the workflow as a whole.
 
