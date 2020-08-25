################################################
Dockerizing Your WDL When You Are Ready For JAWS
################################################

This page will describe how to configure your WDL to use a **docker** image that you have already created.  

Docker images and Conda environments can both be used in WDLs; however, docker is more portable. 

*******************************************
Constraints when Creating your WDL for JAWS
*******************************************

1) one docker image per task - this is a general constraint that `Cromwell <https://cromwell.readthedocs.io/en/stable/>`_ has. 

2) docker image must be registered with `docker hub <https://hub.docker.com>`_  - this is how we have set up the docker backend configuration.


*************************************
Running Your Commands in a Container
*************************************
singularity, docker and shifter are supported depending on which site you are using for JAWS jobs.  

The currently available resources supported by JAWS is

+-----------------------+-------------+
|Available Compute Site | What it Uses|
+=======================+=============+
|JGI (LBLs lawrencium)  | singularity |
+-----------------------+-------------+
|CORI (nersc)           | shifter     |
+-----------------------+-------------+

No matter which site you are using, the WDL is identical.
