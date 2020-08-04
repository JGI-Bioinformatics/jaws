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

   (a) a README.md file describing the workflow
   (b) the WDL 


.. note::
   Keep in mind your README.md will be public and provide users with info on what to expect from the workflow. For ideas on some things to include, see this template :doc:`README.md </Tutorials/suggested_readme_template>`


Now add your WDL to the catalog

.. code-block:: bash

	jaws wdl add <wdl name> <my.wdl> <version> <readme.md>


Use this command to see if your WDL was registered with JAWS 

.. code-block:: bash

   jaws wdl list

