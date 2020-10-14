======================================================
Using Docker Images to Define your Running Environment
======================================================

.. role:: bash(code)
   :language: bash

In this tutorial, I will describe one way docker images can be created and used in your WDL. If you are unfamiliar with docker, please see `docker tutorial <https://scotch.io/tutorials/getting-started-with-docker>`_ or search for the many YouTube tutorials.

.. note::
	As a pre-requisite, you will need a computer with docker installed (Docker Engine - Community).  Installation instructions can be found at `docs.docker.com/install <https://docs.docker.com/install/>`_ or if you have conda installed :bash:`conda install -c conda-forge docker-py`.


Here are the steps we're going to take for this tutorial:

   1. first create an environment with conda - good for testing
   2. make a docker image from conda commands - to be used in the WDL
   3. test a WDL running with the docker containers added


********************************
Clone the Example Repository
********************************
For this tutorial, I will be using the example code from `jaws-example-wdl <https://code.jgi.doe.gov/advanced-analysis/jaws-tutorial-examples.git>`_.
To follow along, do...

.. code-block:: bash

   git clone https://code.jgi.doe.gov/advanced-analysis/jaws-tutorial-examples.git
   cd jaws-tutorial-example/5min_example
   


1) Create your environment
--------------------------

If you haven't done so in the last section, create the conda environment :bash:`bbtools`

.. code-block:: bash

   # create bbtools env
   conda create --name bbtools

   # Now activate your environment.
   source activate bbtools

Install necessary dependencies into your environment

.. code-block:: bash

   # install dependencies
   conda install -y -c bioconda bbmap==38.49
   conda install -y -c bioconda samtools==1.9

3) Testing your environment
---------------------------

For our tutorial, we have a sample wrapper that will become a "task" in the WDL, script.sh. 

.. note :: 
   Each wrapper should write its output to the **current working directory**. You can copy files to other directories after the pipeline has finished.


Try running the script to test environment. You need to be in this directory <your_repo_clone>/jaws-example-wdl/create_wdl_tutorial/. Also, make sure you activated the bbtools environment.

.. code-block:: bash
   
   ./script.sh ../data/sample.fastq.bz2 ../data/sample.fasta

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


When you are convinced the docker image is good, you can register it with `hub.docker.com <hub.docker.com>`_  (you need to make an account first).  When you run a WDL in JAWS, the docker images will be pulled from hub.docker.com. 

.. code-block:: bash

   docker login
   docker push <your_docker_hub_user_name>/bbtools:1.0.0


5) Test your image on cori
--------------------------

Test the docker container on cori.NERSC.gov. You'll need to use the shifter command instead of docker to run your workflow, but the image is the same. More about `shifter at NERSC <https://docs.nersc.gov/development/shifter/how-to-use/>`_.

example:

.. code-block:: bash

   # pull image from hub.docker.com
   shifterimg pull jfroula/bbtools:1.0.0

   # run your wrapper script. notice we are running the script.sh that was saved inside the image
   shifter --image=jfroula/bbtools:1.0.0 script.sh


6) Using the Docker Image in a WDL when Testing
-----------------------------------------------

This subject is handled in more detail on the next page but I will briefly cover it for now. 
Continueing with our :bash:`script.sh` example...

The script.sh that is supplied with the repo has two essential commands: 

.. code-block:: bash
 
   	# align reads to reference contigs
	bbmap.sh in=$READS ref=$REF out=test.sam

	# create a bam file from alignment
	samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam

And it has two inputs :bash:`READS` and :bash:`REF`

.. note:: 
  For testing our WDL, we will inlude the shifter command within the command line section (see below); 
  however, when we are ready to run the WDL in JAWS, the docker image will be removed from the :bash:`command {}` 
  block and added to the :bash:`runtime {}` block, as described in the next section.

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



This sample WDL is also in the repository, called align.wdl.

For a description of what each section of the WDL code does, see the official `WDL docs <https://software.broadinstitute.org/wdl/documentation/quickstart>`_.

