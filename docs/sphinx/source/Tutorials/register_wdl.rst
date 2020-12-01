===================================
How to Register and Share Workflows 
===================================

.. role:: bash(code)
   :language: bash


****************************
Registering your WDL to JAWS
****************************
By adding your workflow to the JAWS catalog, you are making it available to all JAWS users.  This also archives workflows so people can point to older versions for data reproducibility reasons.

There are two required files

   1. README.md file describing the workflow
   2. WDL 


The README.md
-------------

The purpose of the README.md will be to allow the public to run and understand your workflow and even take it to modify for their own purposes. 

The README.md should include things like:

    1. **what does workflow accomplish.**

    2. **description of important tasks/stages.**
       The first two points need to help a potential user understand if they should use this workflow.

    3. **how to recreate running environment.**
       We want to encourage the re-use of code so if docker images were used in the WDL, the location of the dockerfile should be provided (e.g. in hub.docker.com) or a Dockerfile saved inside the developers git repository. Otherwise, explain how the running environment can be re-created so a potential developer could take this WDL and modify it. 
    4. **description of the inputs.**
    5. **how to interpret the outputs.**
    6. **Owner of the workflow.** 
       Who will be responsible for updating this workflow and what is their email.   



Add your WDL to the catalog
---------------------------

.. code-block:: text

    jaws wdl add <wdl name> <my.wdl> <version> <readme.md>


Use this command to see if your WDL was registered with JAWS 

.. code-block:: text

   jaws wdl list

