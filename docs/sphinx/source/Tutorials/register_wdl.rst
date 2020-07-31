===================================
How to Register and Share Workflows 
===================================

.. role:: bash(code)
   :language: bash


****************************
Registering your WDL to JAWS
****************************
By adding your workflow to the JAWS catalog, you are making it available to all JAWS users.  You are also 
guaranteeing a snapshot version that anyone can point to for reproducibility.

.. code-block:: bash

   git clone git@bitbucket.org:berkeleylab/jgi-workflows.git
   cd jgi-workflows
   
   # Create a folder corresponding to the name of your workflow. This will be the public workflow name.
   mkdir <my_public_workflow_name>

Add at least the following two files:

   (a) a README.md file describing the workflow
   (b) the WDL/CWL workflow file named by it's version.  (e.g. :bash:`v2.1.9.wdl`).  You may have multiple WDL files, each corresponding to a different version.  You may also add release notes for each version by adding \*.md files with the same basename as the corresponding WDL file.


If you have two versions of your WDL, you might see something like this in your directory

.. code-block:: bash
   
   2.1.0.md  2.1.0.wdl  2.1.9.wdl	README.md


.. note::
   When you create your own README.md, keep in mind this will be public and provide users with info on what to expect from the workflow. For ideas on some things to include, see this template :doc:`README.md </Tutorials/suggested_readme_template>`

   **(not quite implemented yet)**
   In the :doc:`Catalog </Specifications/catalog_wdls>` of WDLs, you will be able to see the README for each workflow, just click on a name. 


Now push your files to the main repo

.. code-block:: bash

   git add .
   git commit -m 'added my_public_workflow_name'
   git push


Use this command to see if your wdl was registered with JAWS (may take up to 1 hr. to be added to registry)

.. code-block:: bash

   wf list

