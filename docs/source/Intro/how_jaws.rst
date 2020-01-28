==================
How JAWS Works
==================

Pipelines are managed through a workflow engine `(Cromwell) <https://cromwell.readthedocs.io/en/stable/>`_ that requires the tasks in a pipeline to be wrapped in WDL (Workflow Description Language).  WDL provides an intuitive human-readable syntax.

To make running pipelines easier, we have build a robust framework using existing, proven technology that conveniently masks the details of how and where the jobs are run.

#########################
Diagram of how JAWS works
#########################

.. figure:: /Figures/jaws_diagram.png
   :scale: 100%

JAWS is composed of three main parts. 
  
   1) A command line interface (jaws cli)
   2) a workflow engine (Cromwell)
   3) job submission to worker pools (JTM)


JAWS-cli
-----------------
The main purpose of JAWS-cli
  * The JAWS-cli is a commandline interface for the user and gives you a set of commands to submit & monitor jobs, and access logs.
  * JAWS-cli gets everything prepared for Cromwell to take over. For example, it points Cromwell to the wdl, input.json, reference databases, etc.
  * JAWS saves information about runs, like metadata, which version of the pipeline was run, how long it took, etc..; the JAWS-cli gives users access to these records.

Cromwell
----------
The main purpose of Cromwell
  * The Cromwell workflow execution engine that can run WDL workflows locally, in the cloud, or in our case, we've defined a custom backend "JTM"
  * Cromwell parses the WDL, generates individual jobs and dispatches them for execution via the JTM backend.
  * Since Cromwell knows which tasks are dependend on another and which tasks can be run immediately, it can submit tasks to the backend in an efficient manner, parallelizing when possible.

JTM
---
The main purpose of JTM
  * Receive and execute Cromwell tasks from JAWS
  * Share and reuse pools of workers to process tasks from multiple pipelines (reuseing pools means you don't have to wait in the slurm queue again)
  * Hiding the cluster/cloud types and dependencies. The worker(s) can be run on any HPC/Cloud. The backend is maintained by JTM, so a user doesn’t need to consider the backend configurations).
  * Each task in a pipeline can be run on a machine that matches the task's requirements.  For example, an assembly task that needs a lot of memory would be delegated to a large mem 500G machine while another task could be delegated to a lower memory machine.


##################
Technologies used:
##################
- **Authentication:** Globus OAuth
- **Cromwell:** processes workflows described in either WDL `Workflow Description Language <https://software.broadinstitute.org/wdl>`_ or `CWL(Common Workflow Language) <https://www.commonwl.org>`_.
- **Docker, Shifter, Singularity, or conda:** defines run environment
- **JGI Task Manager (JTM):** jobs are relayed to multiple compute clusters; e.g. SLURM-managed Cori (NERSC) & Lawrencium (LBNL) clusters, AWS
- **Globus:** File transfer to/from multiple end-points using GridFTP
- **REST APIs:** Execute workflows either from command line or as a web service 


