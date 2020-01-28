============
Subworkflows
============

.. role:: bash(code)
   :language: bash

.. note::
   This tutorial describes why you would use subworkflows and how to use them in your WDLs. 


A subworkflow is just a WDL that can be called/imported from a main WDL. 

####################
Why use subworkflows
####################
Reasons to use subworkflows
  * re-use code
  * use code that has been expertly vetted

#######
Summary
#######

A subworkflow can be imported into the main WDL with an import statement, for example (:bash:`import "subworkflow.wdl" as subwdl`). It can then be called just like a task would be and it's output can be used by the next task. 

To see the detailed subworkflow description on `cromwell's docs <https://cromwell.readthedocs.io/en/stable/SubWorkflows>`_

Main WDL: main.wdl
------------------

.. figure:: /Figures/main.png

Subworkflow: sub.wdl
---------------------

.. figure:: /Figures/sub.png

##############################
How to Run Subworflows in JAWS
##############################

.. code-block:: bash

   zip -r subworkflows.zip subworkflows.wdl
   jaws -w workflow.wdl --imports subworkflows.zip --inputs input.json

You can try running this example by downloading `subworkflow-helloworld <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/subworkflow-helloworld>`_ repository and running

.. code-block:: bash

   jaws -f main.wdl --imports sub.wdl.zip -i fake.json
