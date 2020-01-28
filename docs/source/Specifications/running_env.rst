###################################################
How to Configure your WDL to use a Docker Container
###################################################

This page will describe how to configure your WDL to use a **docker** image that you have already created.  

Docker images and Conda environments can both be used in WDLs; however, docker is more portable. 

To create a "running environment", you will have to install the libraries and software dependencies required to run your application. To see how to create an environment, you can see the section "Create your environment" at :doc:`Add workflows to JAWS (the whole process) </Tutorials/whole_wdl_process>`.  

*******************************************
Constraints when Creating your WDL for JAWS
*******************************************

1) one docker image per task - this is a general constraint that cromwell has. 

2) docker image must be registered with `docker hub <https://hub.docker.com>`_  - this is how we have set up the docker backend configuration.


*************************************
Running Your Commands in a Container
*************************************
singularity, docker and shifter are supported depending on which cluster you are using for JAWS jobs.  

Example clusters are 

* LBLs lawrencium (singularity)

* AWS (docker)

* NERSC-cori (shifter)

No matter which cluster you are using, the WDL uses similar syntax:

Lawrencium  & AWS
-----------------
For AWS or Lawrencium, you are using singularity or docker, respectively, behind the scenes to run your WDL, so your WDL configuration would look like:

.. code-block:: bash

    command {
	  # all commands are run in the container
	}

	runtime {	
		docker: "jfroula/jaws-blastplus:1.0.16"
	}


Cori
----
For Cori you are using shifter behind the scenes to run your task so your WDL configuration would look like:

.. code-block:: bash

    command {
	  # all commands are run in the container
	}

	runtime {	
		docker: "jfroula/jaws-blastplus:1.0.16"
		docker_img: "jfroula/jaws-blastplus"
		docker_tag: "1.0.16"
	}


.. note:: **Other Options**
	Of course there are ways around these constraints that are not too bad like including the docker or shifter command in the `command` section of the WDL.  Or you can use conda environments instead of docker images.

