==============================
Registering Your WDL with JAWS
==============================

.. role:: bash(code)
   :language: bash

***************
Overall Process
***************

1) **Write a WDL of your pipeline**.
   When you develop your WDL, you can test using any installation of cromwell.jar (and probably not JAWS).

2) **Create a git repository that includes**

   * README.md
   * WDL file
   * inputs json file for testing
   * any test data
   * a configuration file if you need to test tasks that will require sbatch

   Use this repository for testing your WDL.

3) **Clone the official jgi_workflow repository**.
   Once you finish testing the WDL, clone the jgi_workflows repository and add a folder, named something reasonable since this will be the official workflow name.
   Add only the WDL and README.md file.  The prefix of the WDL file will be the version of the workflow.  So for instance, if you created a folder named
   "jgi_metagenome_assembly" and added v1.0.0.wdl to it, the name you would see in JAWS would be :bash:`jgi_metagenome_assembly/v1.0.0`.

.. note::
   You may have to get access to the repository from (<Ed Kirton>eskirton@lbl.gov)

.. code-block:: bash

   # clone jgi_workflows**
   git clone git@gitlab.com:eskirton/jgi_workflows.git
   or
   git clone https://gitlab.com/eskirton/jgi_workflows.git

3) **Register the workflow**
   Push your changes to the jgi_workflows repository to the main JAWS repository https://gitlab.com/eskirton/jgi_workflows.
   This jgi_workflows repo will be seen by JAWS (refreshes every hour or you can run :bash:`jaws --pull`)

.. code-block:: bash

   # add your WDL folder to jgi_workflows*
   git add .
   git commit -m 'adding leo_dap'
   git push

your workflow will now be available to anyone using the system (i.e. :bash:`wf list`) and visible in the JAWS catalog.

.. note::
   You don't need to register your WDL in JAWS to test it on JAWS. You can run JAWS by pointing to your wdl and inputs file :bash:`jaws <wdl> <inputs>`





