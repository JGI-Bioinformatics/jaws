==========================
Creating Docker Containers
==========================

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
To follow along, do:

.. code-block:: text

   git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
   cd jaws-tutorial-examples/5min_example


*******************
Create docker image
*******************

Next we'll describe how to create a Dockerfile and register it with `hub.docker.com <https://docs.docker.com/docker-hub/>`_.  But first create an account and click on "Create a Repository". In the space provided, enter a name for your container, that doesn't have to exist yet, like :bash:`aligner-bbmap`.  You will push a docker image to this name after you create it in the next steps.

To make the Dockerfile, you can use the same commands you used for the conda environment.  Notice that it is good practice to specify the versions when installing software like I have done in the example Dockerfile. Of course, you can drop the versions altogether to get the latest version but the Dockerfile may not work out-of-the-box in the future due to version conflicts.

.. note::
    It is helpful, when creating the Dockerfile to test each command (i.e. apt-get, wget, conda install, etc) manually, inside an empty docker container. Once everything is working, you can copy the commands to a Dockerfile.

This docker command will create an interactive container with an ubuntu base image.  You can start installing stuff as root.

.. code-block:: text

    docker run -it ubuntu:latest /bin/bash


Here is an example Dockerfile (provided in 5min_example). We will create a container from it.

.. code-block:: text

    FROM ubuntu:22.04

    # Install stuff with apt-get
    RUN apt-get update && apt-get install -y wget bzip2 \
    	&& rm -rf /var/lib/apt/lists/*

    # Point to all the future conda installations you are going to do
    ENV CONDAPATH=/usr/local/bin/miniconda3
    ENV PATH=$CONDAPATH/bin:$PATH

    # Install miniconda
    # There is a good reason to install miniconda in a path other than its default.
    # The default intallation directory is /root/miniconda3 but this path will not be
    # accessible by shifter or singularity so we'll install under /usr/local/bin/miniconda3.
    RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh \
        && bash ./Miniconda3*.sh -b -p $CONDAPATH \
        && rm Miniconda3*.sh

    # Install software with conda
    RUN conda install -c bioconda bbmap==38.84 samtools==1.11 \
    	&& conda clean -afy

    # This will give us a workingdir within the container (e.g. a place we can mount data to)
    WORKDIR /bbmap

    # Move script into container.
    # Notes that it is copied to a location in your $PATH
    COPY script.sh /usr/local/bin/script.sh


| **Build the image and upload to hub.docker.com**
| You need to use your docker hub user name to tag the image when you are building it (see below).

.. code-block:: text

   # create a "Build" directory and create docker container from there so its a small image. Its good practice to always create an image in
   # a directory containing only the required files, otherwise the container will also include them and could be very large.
   mkdir build
   cp script.sh Dockerfile build/
   cd build
   docker build --tag <your_docker_hub_user_name>/aligner-bbmap:1.0.0 .
   cd ../


**Test that the example script runs in the docker container**

.. code-block:: text

   # use your image name
   docker run <your_docker_hub_user_name>/aligner-bbmap:1.0.0 script.sh

   # if you are in the root of the 5min_example directory, then try re-running the script with data.
   docker run --volume="$(pwd)/../data:/bbmap" <your_docker_hub_user_name>/aligner-bbmap:1.0.0 script.sh sample.fastq.bz2 sample.fasta

   # Notice script.sh is found because it was set in PATH in the Dockerfile and
   # the two inputs are found because the data directory is mounted to /bbmap (inside container) where the script runs.


When you are convinced the docker image is good, you can register it with `hub.docker.com <hub.docker.com>`_  (remember to make an account first).  When you run a WDL in JAWS, the docker images will be pulled from hub.docker.com.

.. code-block:: text

   docker login
   docker push <your_docker_hub_user_name>/aligner-bbmap:1.0.0

Now your image is available on any site i.e. cori, jgi, tahoma, aws, etc.  And although you can manually pull your image using `shifter pull <https://docs.nersc.gov/development/shifter/how-to-use/#downloading-shifter-images-to-nersc>`_, `singularity pull <https://docs.sylabs.io/guides/3.2/user-guide/cli/singularity_pull.html>`_, or `docker pull <https://docs.docker.com/engine/reference/commandline/pull/>`_, JAWS will do this for you (but cromwell wont).


***********************
Test your image on cori
***********************

Besides your docker-machine, it is useful to test your image on CORI since you will likely be running your WDL there at some point.  There are certain aspects of the docker container that will work on your docker-machine but won't on another site, like cori. This is because shifter or singularity behave differently than docker.

To test the docker container on :bash:`cori.nersc.gov`. You'll need to use the shifter command instead of docker to run your workflow, but the image is the same. More about `shifter at NERSC <https://docs.NERSC.gov/programming/shifter/how-to-use/>`_.

example:

.. code-block:: text

   # pull image from hub.docker.com
   shifterimg pull <your_docker_hub_user_name>/aligner-bbmap:1.0.0

   # clone the repo on cori
   git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
   cd jaws-tutorial-examples/5min_example

   # run your wrapper script. notice we are running the script.sh that was saved inside the image
   shifter --image=<your_docker_hub_user_name>/aligner-bbmap:1.0.0 ./script.sh ../data/sample.fastq.bz2 ../data/sample.fasta


*******
The WDL
*******
The :bash:`script.sh` that is supplied with the repo has two essential commands:

.. code-block:: text

    # align reads to reference contigs
    bbmap.sh in=$READS ref=$REF out=test.sam

    # create a bam file from alignment
    samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam

It would make sense to have both commands inside one task of the WDL because they logically should be run together.  However, for an excersise, we will have the two commands become two tasks.  The output from the first command is used in the second command, so in our WDL example, we can see how tasks pass information.

See an example of the finished WDL :bash:`align_final.wdl` and its :bash:`input.json`` file

.. raw:: html

   <details>
   <summary style="color: #448ecf";>align_final.sh</summary>

.. code-block:: text

    version 1.0

    workflow bbtools {
        input {
            File reads
            File ref
        }

        call alignment {
           input: fastq=reads,
                  fasta=ref
        }
        call samtools {
           input: sam=alignment.sam
       }
    }

    task alignment {
        input {
            File fastq
            File fasta
        }

        command {
            bbmap.sh in=~{fastq} ref=~{fasta} out=test.sam
        }

        runtime {
            docker: "jfroula/aligner-bbmap:2.0.2"
            time: "00:10:00"
            memory: "5G"
            cpu: 1
        }

        output {
           File sam = "test.sam"
        }
    }

    task samtools {
        input {
            File sam
        }

        command {
           samtools view -b -F0x4 ~{sam} | samtools sort - > test.sorted.bam
        }

        runtime {
            docker: "jfroula/aligner-bbmap:2.0.2"
            time: "00:10:00"
            memory: "5G"
            cpu: 1
        }

        output {
           File bam = "test.sorted.bam"
        }
    }

.. raw:: html

    </details>
    <br>

.. raw:: html

   <details>
   <summary style="color: #448ecf";>inputs.json</summary>

.. code-block:: text

    {
        "bbtools.reads": "../data/sample.fastq.bz2",
        "bbtools.ref": "../data/sample.fasta"
    }

.. raw:: html

    </details>
    <br><br>

.. note::
    Singularity, docker, or shifter can be prepended to each command for testing (see align_with_shifter.sh); however,
    this wouldn't be appropriate for a finished "JAWSified" WDL because you loose portability.  The final WDL should have the docker image name put inside the :bash:`runtime {}` section.

This may be helpful when testing & debugging so I've included an example where shifter is prepended to each command.

.. raw:: html

   <details>
   <summary style="color: #448ecf";>align_with_shifter.sh</summary>

.. code-block:: text

    version 1.0

    workflow bbtools {
        input {
            File reads
            File ref
        }

        call alignment {
           input: fastq=reads,
                  fasta=ref
        }
        call samtools {
           input: sam=alignment.sam
       }
    }

    task alignment {
        input {
            File fastq
            File fasta
        }

        command {
            shifter --image=jfroula/aligner-bbmap:2.0.1 bbmap.sh in=~{fastq} ref=~{fasta} out=test.sam
        }

        output {
           File sam = "test.sam"
        }
    }

    task samtools {
        input {
            File sam
        }

        command {
           shifter --image=jfroula/aligner-bbmap:2.0.1 samtools view -b -F0x4 ~{sam} | shifter --image=jfroula/aligner-bbmap:2.0.1 samtools sort - > test.sorted.bam
        }

        output {
           File bam = "test.sorted.bam"
        }
    }

You would run this WDL on Cori with the following command.

.. code-block:: text

    java -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run align_with_shifter.wdl -i inputs.json

.. raw:: html

   </details>
   <br><br>


***************************************************
The Docker Image Should be in the runtime{} Section
***************************************************

Everything in the :bash:`command{}` section of the WDL will run inside a docker container if you've added docker to the :bash:`runtime{}` section.
Now your WDL has the potential to run on a machine with shifter, singularity, or docker.  JAWS will take your docker image and run it appropriately as singularity, docker or shifter.  If you run the WDL with the cromwell command on a shifter or singularity machine, you need to supply a :bash:`cromwell.conf` file, explained shortly.

See :bash:`align_final.wdl`:

.. code-block:: text

    runtime {
        docker: "jfroula/aligner-bbmap:1.2.1"
    }

.. _run with conf:

*******************************
Run the Final WDL with Cromwell
*******************************

On a docker machine
-------------------

You can now run the final WDL:

.. code-block:: text

    conda activate bbtools  # you need this for the cromwell command only
    cromwell run align_final.wdl -i inputs.json


On Cori
-------
You'll have to include a cromwell.conf file in the command because it is the config file that knows whether to run the image, supplied in the :bash:`runtime{}` section, with docker, singularity, or shifter.  You don't need to supply a cromwell.conf file in the above cromwell command because docker is default.

The cromwell.conf file is used to:

1. override cromwell's default settings
2. tells cromwell how to interpret the WDL (i.e. use shifter, singularity, etc)
3. specifies the backend to use (i.e. local, slurm, aws, condor, etc)

.. note::

	JAWS takes care of the cromwell.conf for you.


Here you can find the config files: `jaws-tutorials-examples/config_files <https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples/-/tree/master/config_files>`_.


.. code-block:: text

    java -Dconfig.file=../config_files/cromwell_cori.conf \
         -Dbackend.providers.Local.config.dockerRoot=$(pwd)/cromwell-executions \
         -Dbackend.default=Local \
         -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run align_final.wdl -i inputs.json

where

    :bash:`-Dconfig.file`
    points to a cromwell conf file that is used to overwrite the default configurations

    :bash:`-Dbackend.providers.Local.config.dockerRoot`
    this overwrites a variable 'dockerRoot' that is in cromwell_cori.conf so that cromwell will use your own current working directory to place its output.

    :bash:`-Dbackend=[Local|Slurm]`
    this will allow you to choose between the Local and Slurm backends. With slurm, each task will have it's own sbatch command (and thus wait in queue).

*********************************
Understanding the Cromwell Output
*********************************
Cromwell output is:

1. files created by the workflow
2. the stdout/stderr printed to screen

**1. Where to find the output files**

Cromwell saves the results under a directory called :bash:`cromwell-executions`.  And under here, there is a unique folder name representing one WDL run.

.. figure:: /Figures/crom-exec.svg
    :scale: 100%

Each task of your workflow gets run inside the :bash:`execution` directory so it is here that you can find any output files including the stderr, stdout & script file.

Explaination of cromwell generated files

.. raw:: html

   <details>
   <summary style="color: #448ecf";>stderr</summary>
   <p class="textborder">
   The stderr from any of the commands/scripts in your task should be in this file.
   </details>

   <details>
   <summary style="color: #448ecf";>stdout</summary>
   <p class="textborder">
   The stdout from all the commands/scripts in your task should be in this file. Not all scripts send errors to stderr as they should so you will find them in here instead.
   </details>

   <details>
   <summary style="color: #448ecf";>script</summary>
   <p class="textborder">
   The script file is run by the script.submit file.  It contains all the commands that you supplied in the commands{} section of the WDL, as well as cromwell generated code that creates the stderr, stdout, and rc files.
   </details>

   <details>
   <summary style="color: #448ecf";>script.submit</summary>
   <p class="textborder">
   This file contains the actual command that cromwell ran.  If the file was created by JAWS, there is one more step before "script" gets run.
   <br><br>script.submit -> dockerScript -> script
   </details>

   <details>
   <summary style="color: #448ecf";>rc</summary>
   <p class="textborder">
   This file contains the return code for the commands{} section of the WDL.  One thing to remember is that the return code used for the rc file is from your last command run.  And so if a command fails but the last command succeeded, the return code would be 0, unless you used "set -e" which forces an exit upon the first error.
   </details>
   <br>

These files are only seen in JAWS

.. raw:: html

   <details>
   <summary style="color: #448ecf";>stdout.submit</summary>
   <p class="textborder">
   This file is created by script.submit and not by the script file and the content is not useful for debugging your task.
   </details>

   <details>
   <summary style="color: #448ecf";>stderr.submit</summary>
   <p class="textborder">
   This file is created by script.submit and not by the script file which means there may be some useful error messages.  If there was a problem upstream of the task even starting, the error should be in this file.
   </details>

   <details>
   <summary style="color: #448ecf";>dockerScript</summary>
   <p class="textborder">
   This file is created by script.submit and runs the script file.
   <br><br>script.submit -> dockerScript -> script
   </details>
   <br>

**2. Cromwell's stdout**

When you ran :bash:`align_with_shifter.wdl` with cromwell above, observe these lines in the output.

1. the bash bbmap.sh and samtools commands that were run
2. paths to the output files from the workflow
3. you should see WorkflowSucceededState
4. copy a path from one of the output execution directories. Notice the cromwell generated files and your :bash:`.sam` or :bash:`.bam` output is there.
5. :bash:`Call-to-Backend` shows that we are running on local backend (default)

.. note::
  You won't have access to this same cromwell standard output when you run through JAWS.  The same information can be found in different ways.


Limitations when using docker
-----------------------------
1. One docker image per task - this is a general constraint that Cromwell has.
2. The docker image must be registered with docker hub - this is how we have set up the docker backend configuration.
3. A `sha256` tag must be used instead of some custom tag (i.e v1.0.1) for call-caching to work.

    To find the `sha256` tag, you can use:

    .. code-block:: text

        # on a docker-machine
        docker images --digests | grep <your_docker_hub_user_name>

        # on a shifter-machine
        shifterimg lookup ubuntu:16.04

    The version tag (16.04) can be replaced by the sha256 tag.

    .. code-block:: text

        runtime {
            docker: "ubuntu@sha256:20858ebbc96215d6c3c574f781133ebffdc7c18d98af4f294cc4c04871a6fe61"
        }

    You can interactively go into a container from shifter by

    .. code-block:: text

        shifter --image=id:20858ebbc96215d6c3c574f781133ebffdc7c18d98af4f294cc4c04871a6fe61
        or
        shifter --image=ubuntu:16.04

