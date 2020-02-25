==================
How JAWS Works
==================

Pipelines are managed through a workflow engine `(Cromwell) <https://Cromwell.readthedocs.io/en/stable/>`_ that requires the tasks in a pipeline to be wrapped in WDL (Workflow Description Language).  WDL provides an intuitive human-readable syntax.

To make running pipelines easier, we have build a robust framework using existing, proven technology that conveniently masks the details of how and where the jobs are run.


JAWS is composed of three main parts. 
  
   1l a command line interface (JAWS Client)
   2) a workflow engine (Cromwell)
   3) job submission to worker pools (JTM)

#########################
Diagram of how JAWS works
#########################
The main takeaways here are 

  * all commands are from the command line ("client")
  * the "central" installation of the JAWS server coordinates with sites (e.g. LBNL or NERSC) where the pipelines are run. 
  * GLOBUS transfers all your files from your data source to the computing-site where Cromwell will actually run. 
  * Cromwell will run the pipeline at the computing-site (e.g. LBNL or NERSC).
  * the JTM (JGI Task Manager) serves as the backend to Cromwell and handles the running of the jobs. 


.. figure:: /Figures/JAWS-Arch.svg
   :scale: 100%



###################################################
Details of JAWS Implementation at LBL IT and NERSC  
###################################################
This diagram emphasizes the differences between the two site implementations and how everything communicates.  A few redundant features were excluded due to space limitations.


.. figure:: /Figures/JAWS-System.svg
   :scale: 100%

JAWS-Client
-----------------
The main purpose of JAWS-cli
  * the JAWS-Client is a command-line interface for the user and gives you a set of commands to submit & monitor jobs, and access logs.
  * the JAWS-Client gets everything prepared for Cromwell to run the workflow. For example, it points Cromwell to the WDL, input JSON, reference databases, etc.
  * JAWS saves information about runs, like which version of the pipeline was run, how long it took, etc..; the CLI provides a way for the user to retrieve these records

Cromwell
----------
The main purpose of Cromwell
  * the Cromwell workflow execution engine that can run WDL workflows locally, in the cloud, or in our case, we've defined a custom backend "JTM"
  * Cromwell parses the WDL, generates individual jobs and dispatches them for execution via the JTM backend.
  * since Cromwell knows which tasks are dependend on another and which tasks can be run immediately, it can submit tasks to the backend in an efficient manner, parallelizing when possible.

JTM
---
The main purpose of JTM
  * receive tasks from Cromwell and execute them on requestable computing resources (e.g. cluster)
  * share and reuse pools of workers to process tasks from multiple pipelines (reuseing pools means you don't have to wait in the cluster scheduler's queue again)
  * hiding the cluster/cloud types and dependencies. The worker(s) can be run on any HPC/Cloud. The backend is maintained by JTM, so a user doesn't need to consider the backend configuration.  No changes to the WDL are required to run it on a different compute-cluster.
  * each task in a pipeline can be run on a machine that matches the task's requirements.  For example, an assembly task that needs a lot of memory would be delegated to a large mem 500G machine while another task could be delegated to a lower memory machine.


##################
Technologies used:
##################
- **Authentication:** Globus OAuth
- **Cromwell:** processes workflows described in either WDL `Workflow Description Language <https://software.broadinstitute.org/WDL>`_ or `CWL(Common Workflow Language) <https://www.commonwl.org>`_.
- **Docker, Shifter, Singularity, or conda:** defines run environment
- **JGI Task Manager (JTM):** jobs are relayed to multiple compute clusters; e.g. SLURM-managed Cori and Lawrencium clusters, AWS
- **Globus:** File transfer to/from multiple end-points using GridFTP
- **REST APIs:** Execute workflows either from command line or as a web service 


