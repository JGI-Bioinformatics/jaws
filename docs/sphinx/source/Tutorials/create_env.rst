======================================================
Create a Running Environment Using a Docker Container
======================================================

.. role:: bash(code)
   :language: bash

*******
Summary
*******

In this tutorial, I will describe one way docker images can be created and used in your WDL. If you are unfamiliar with docker, please see `docker tutorial <https://scotch.io/tutorials/getting-started-with-docker>`_ or search for the many YouTube tutorials.

.. note::
    As a pre-requisite, you will need a computer with docker installed (Docker Engine - Community).  Installation instructions can be found at `docs.docker.com/install <https://docs.docker.com/install/>`_ or if you have conda installed :bash:`conda install -c conda-forge docker-py`.  


Here are the steps we're going to take for this tutorial:
   1. make a docker image from the same commands you used for the conda environment (last tutorial)
   2. run a WDL that is using your docker container


****************************
Clone the Example Repository
****************************
For this tutorial, I will be using the example code from `jaws-tutorial-examples <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git>`_.
To follow along, do...

.. code-block:: text

   git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
   cd jaws-tutorial-example/5min_example
   

*******************
Create docker image
*******************

Next we'll describe how to create a Dockerfile and register it with `hub.docker.com <https://docs.docker.com/docker-hub/>`_.

.. note::
    you'll have to create an account and an empty repository with hub.docker.com first.

To make the Dockerfile, you can use the same commands you used for the conda environment.  Notice that it is good practice to specify the versions for each software like I have done in the example Dockerfile. Of course, you can drop the versions altogether to get the latest version but the Dockerfile may not work out-of-the-box in the future due to version conflicts.


The Dockerfile (provided in example) looks like 

.. code-block:: text

    FROM ubuntu:16.04

    # install stuff with apt-get
    RUN apt-get update && apt-get install -y wget bzip2
    
    # Install miniconda
    # There is a good reason to install miniconda in a path other than its default.
    # The default intallation directory is /root/miniconda3 but this path will not be
    # accessible by shifter or singularity so we'll install under /usr/local/bin/miniconda3.
    RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh \
        && bash ./Miniconda3-4.5.11-Linux-x86_64.sh -b -p /usr/local/bin/miniconda3 \
        && rm Miniconda3-4.5.11-Linux-x86_64.sh
    
    # point to all the future conda installations you are going to do
    ENV PATH=/usr/local/bin/miniconda3/bin:$PATH
    
    # install software
    RUN conda install -c bioconda bbmap==38.84
    RUN conda install -c bioconda samtools==1.11
    
    # this will give us a workingdir within the container (e.g. a place we can mount data to) 
    WORKDIR /bbmap
    
    # move script into container
    COPY script.sh /usr/local/bin/script.sh



| **Build the image and upload to hub.docker.com**
| You need to use your docker hub user name to tag the image when you are building it (see below).

.. code-block:: text

   # create a "Build" directory and create docker container from there so its a small image. Its good practice to always create an image in 
   # a directory containing only the required files, otherwise the container will also include them and could be very large.
   mkdir build 
   cp script.sh Dockerfile build/
   cd build
   docker build --tag <your_docker_hub_user_name>/bbtools:1.0.0 .
   cd ../


**Test that the script runs in the docker container**

.. code-block:: text

   # use your image name
   docker run jfroula/bbtools:1.0.0 script.sh
 
   # if you are in the root of the 5min_example directory, then try re-running the script with data.
   docker run --volume="$(pwd)/../data:/bbmap" jfroula/bbtools:1.0.0 script.sh sample.fastq.bz2 sample.fasta

   # Notice script.sh is found because it was set in PATH in the Dockerfile and
   # the two inputs are found because the data directory is mounted to /bbmap (inside container) where the script runs.



When you are convinced the docker image is good, you can register it with `hub.docker.com <hub.docker.com>`_  (remember to make an account first).  When you run a WDL in JAWS, the docker images will be pulled from hub.docker.com. 

.. code-block:: text

   docker login
   docker push <your_docker_hub_user_name>/bbtools:1.0.0


***********************
Test your image on cori
***********************

Test the docker container on cori.NERSC.gov. You'll need to use the shifter command instead of docker to run your workflow, but the image is the same. More about `shifter at NERSC <https://docs.NERSC.gov/programming/shifter/how-to-use/>`_.

example:

.. code-block:: text

   # pull image from hub.docker.com
   shifterimg pull jfroula/bbtools:1.0.0

   # clone the repo on cori
   git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
   cd jaws-tutorial-example/5min_example

   # run your wrapper script. notice we are running the script.sh that was saved inside the image
   shifter --image=<your_docker_hub_user_name>/bbtools:1.0.0 ./script.sh ../data/sample.fastq.bz2 ../data/sample.fasta


*******************
Add Docker to a WDL
*******************
The :bash:`script.sh` that is supplied with the repo has two essential commands: 

.. code-block:: text
 
    # align reads to reference contigs
    bbmap.sh in=$READS ref=$REF out=test.sam

    # create a bam file from alignment
    samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam

It would make sense to have both commands inside one task of the WDL because they logically should be run together.  However, for an excersise, we will have the two commands become two tasks.  The output from the first command is used in the second command, so in our WDL example, we can see how tasks pass information.

The docker command (or shifter if you are on cori) can be appended to each command for testing. This wouldn't be appropriate for a finished "JAWSified" WDL because you loose portability.  The final WDL should have the docker image name put inside the :bash:`runtime {}` section.


See file align_with_shifter.sh

.. code-block:: text

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



*****************************
Running the WDL Using Shifter
*****************************

| Now when you run align_with_shifter.wdl, you don't need your conda environment.
| (this will only work on cori which supports shifter)

.. code-block:: text

    java -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run align_with_shifter.wdl -i inputs.json


**********************************************
Move the Docker Image to the runtime{} Section
**********************************************

.. note::
    To get a description of the runtime section, see :ref:`requesting-workers`.

After shifter is removed from the :bash:`command{}` block, add :bash:`docker:` inside the :bash:`runtime{}` block to each of the tasks in the WDL. Now, all the code inside :bash:`commands{}` will be run inside a container.

See align.wdl

.. code-block:: text

    runtime {
        docker: "jfroula/bbtools:1.2.1"
    }

.. _run with conf:

*************************************
Run with Docker Inside the runtime{}
*************************************

To run again you have to use a slightly different command which overwrites the default :bash:`dockerRoot` path so it points to your current working directory. 

This also has to be run on **cori** since the config file uses shifter to run the container. This could instead be configured with the docker command if you wanted to test on your laptop.  Here's an example of a `docker configuration <https://github.com/broadinstitute/cromwell/blob/develop/cromwell.example.backends/Docker.conf>`_


.. code-block:: text

    java -Dconfig.file=cromwell_cori.conf \
         -Dbackend.providers.Local.config.dockerRoot=$(pwd)/cromwell-executions \
         -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run jgi_meta.jaws.wdl -i inputs.json

where 

    :bash:`-Dconfig.file` 
    points to a cromwell conf file that is used to overwrite the default configurations

    :bash:`-Dbackend.providers.Local.config.dockerRoot`
    this overwrites a variable 'dockerRoot' that is in cromwell_cori.conf so that cromwell will use your own current working directory to place its output.


