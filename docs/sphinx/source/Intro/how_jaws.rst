====================
What is JAWS
====================

.. role:: bash(code)
  :language: bash

`(see video about the basics of JAWS) <https://youtu.be/85lJFvGFVpE>`_

JAWS is a multi-site workflow manager that uses the `Cromwell <https://Cromwell.readthedocs.io/en/stable/>`_ workflow engine. Some main directives of JAWS are to make running of bioinformatics workflows easier, foster collaboration between users of the system, and make it possible to move workloads across different resources.

JAWS is composed of four main parts:
  
	1) a command line interface: "Jaws Client"
	2) a centralized orchestration service: "Jaws Central", administering runs to multiple sites.
	3) a site service that wraps the workflow engine, like Cromwell, and is installed on a compute site.
	4) a job submission manager which submits jobs to worker pools: "JTM" and is designed to efficiently use the high performance computing (HPC) resources.

#########################
Diagram of how JAWS works
#########################
Below is a diagram of the JAWS architecture. Note that there is some duplication of processes that is meant to demonstrate that "site" can be installed at multiple sites.   

The main takeaways here are:

  * All the commands are from the command line and handled by :bash:`Jaws Client`.
  * The :bash:`Jaws Central` is a server that coordinates which compute-site (e.g. LabIT or NERSC) the pipeline is run. 
  * `GLOBUS <https://globus.org/>`_ transfers all your files from your data source to the computing-site where Cromwell will actually run. 
  * Cromwell is the workflow engine that will run the pipeline at the compute-site.
  * The JTM (JGI Task Manager) (:bash:`JTM`) serves as the backend to Cromwell and handles the running of the jobs on a HPC cluster. 

.. figure:: /Figures/jaws_architecture-Architecture_v2_8_5.svg
   :scale: 100%

Click on the image to enlarge   

#################################################
Details of JAWS Implementation at LBL IT and CORI  
#################################################
This diagram emphasizes the differences between the two site implementations and how everything communicates. A few redundant features were excluded due to space limitations.

.. figure:: /Figures/JAWS-System.svg
   :scale: 100%

Click on the image to enlarge  

JAWS Overall Workflow Processing
--------------------------------
The user interfaces only with the :bash:`jaws-client`. The :bash:`jaws-client` communicates with :bash:`jaws-central` to move data to the target site and hands over the workflow executions to the respective :bash:`jaws-site` service which in turn runs the workflow to completion and relays the status back to :bash:`jaws-central`. Globus is used as a transfer mechanism between a central data storage location and target sites. The execution of workflows by :bash:`jaws-site` is orchestrated by Cromwell. 

jaws-client
-----------
:bash:`jaws-client` is a command-line interface for the user and interacts with the central service using defined APIs. :bash:`jaws-client` offers commands to submit and monitor workflows. :bash:`jaws-central` saves metadata about runs, for example, which version of the pipeline was run, runtime statistics, which datasets were processed, etc.

Cromwell
----------
`Cromwell <https://Cromwell.readthedocs.io/en/stable/>`_ is responsible for executing the commands in a workflow. The tasks are executed on a user-defined backend, i.e. :bash:`jaws-jtm`.

JTM (jaws-jtm)
--------------
The main purpose of the JAWS JTM (JGI Task Manager) is to receive tasks from Cromwell and execute them on a compute resource (e.g. HPC cluster). Cromwell sends the workflow tasks to the workers running on the HPC cluster via JTM. JTM accomplishes this by using RabbitMQ message broker.  It acts as an abstraction layer between :bash:`jaws-site` and different resources (different clusters, eventually cloud resources).

Globus
------
`GLOBUS <https://globus.org/>`_ transfers all your files from your data source to the computing-site where Cromwell runs.

##################
Technologies used:
##################
- **Authentication:** Globus OAuth
- **Cromwell:** processes workflows described in either WDL `Workflow Description Language <https://software.broadinstitute.org/WDL>`_.
- **Docker, Shifter, Singularity** defines run environment
- **JGI Task Manager (JTM):** jobs are relayed to multiple compute clusters; for example, Cori & LBNL
- **Globus:** File transfer to/from multiple end-points using GridFTP
- **REST APIs:** multiple JAWS components communicate by REST
- **RabbitMQ:** Message broker used to communicate workflow tasks between Cromwell and the JTM workers running on the compute cluster.


