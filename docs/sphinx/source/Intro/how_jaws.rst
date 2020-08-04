====================
What is JAWS
====================

.. role:: bash(code)
  :language: bash

`(see video about the basics of JAWS) <https://youtu.be/85lJFvGFVpE>`_

JAWS is a multi-site workflow manager, using the `Cromwell <https://Cromwell.readthedocs.io/en/stable/>`_ workflow engine. The main purpose of JAWS is to make running of bioinformatics workflows easier, to foster collaboration between users of the system, and to make it possible to move workloads across different resources.

JAWS is composed of four main parts:
  
	1) a command line interface (:bash:`jaws-client`) 
	2) a centralized orchestration service (:bash:`jaws-central`), managing N sites.
	3) a site service wrapping the workflow engine (:bash:`jaws-site`)
	4) job submission to worker pools, to make workloads amenable to HPC resources (:bash:`jaws-jtm`)



#########################
Diagram of how JAWS works
#########################
The main takeaways here are 

  * All the commands are from the command line and handled by (:bash:`jaws-client`).
  * The :bash:`jaws-central` is a server that coordinates which site (e.g. LBNL or NERSC) the pipeline is run. 
  * GLOBUS transfers all your files from your data source to the computing-site where Cromwell will actually run. 
  * Cromwell will run the pipeline at the computing-site.
  * The JTM (JGI Task Manager) (:bash:`jaws-jtm`) serves as the backend to Cromwell and handles the running of the jobs. 


.. figure:: /Figures/JAWS-Arch.svg
   :scale: 100%

|
|


###################################################
Details of JAWS Implementation at LBL IT and NERSC  
###################################################
This diagram emphasizes the differences between the two site implementations and how everything communicates.  A few redundant features were excluded due to space limitations.


.. figure:: /Figures/JAWS-System.svg
   :scale: 100%


JAWS Overall Workflow Processing
--------------------------------
The user interfaces only with the :bash:`jaws-client`. The :bash:`jaws-client` communicates with :bash:`jaws-central` to move data to the target site and hands over the workflow executions to the respective :bash:`jaws-site` service which in turn runs the workflow to completion and relays the status back to :bash:`jaws-central`. Globus is used as a transfer mechanism between a central data storage location and target sites. The execution of workflows by :bash:`jaws-site` is orchestrated by Cromwell, which is an external dependency.


jaws-client
-----------
:bash:`jaws-client` is a command-line interface for the user and interacts with the central service using defined APIs. :bash:`jaws-client` offers commands to submit and monitor workflows. :bash:`jaws-central` saves metadata about runs, for example, which version of the pipeline was run, runtime statistics, which datasets were processed, etc

Cromwell
----------
Cromwell is responsible for executing the commands in a workflow. The tasks are executed on a user-defined backend, i.e. :bash:`jaws-jtm`.

JTM (jaws-jtm)
--------------
The main purpose of the JAWS JTM (JGI Task Manager) is to receive tasks from Cromwell and execute them on available computing resources (e.g. HPC cluster). Cromwell sends the workflow tasks to the workers running on the HPC cluster via JTM. JTM accomplishes this by using RabbitMQ message broker.  It acts as an abstraction layer between :bash:`jaws-site` and different resources (different clusters, eventually cloud-like resources).

Globus
------
GLOBUS transfers all your files from your data source to the computing-site where Cromwell will actually run.

##################
Technologies used:
##################
- **Authentication:** Globus OAuth
- **Cromwell:** processes workflows described in either WDL `Workflow Description Language <https://software.broadinstitute.org/WDL>`_.
- **Docker, Shifter, Singularity** defines run environment
- **JGI Task Manager (JTM):** jobs are relayed to multiple compute clusters; for example, Cori 
- **Globus:** File transfer to/from multiple end-points using GridFTP
- **REST APIs:** multiple JAWS components communicate by REST
- **RabbitMQ:** Message broker used to communicate workflow tasks between Cromwell and the JTM workers running on the compute cluster.


