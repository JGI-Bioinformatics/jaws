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

A subworkflow can be imported into the main WDL and then treated as a task. 

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

You can try running this example by downloading `subworkflow-helloworld <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/subworkflow-helloworld>`_ as a tar or zip or you can clone :bash:`https://code.jgi.doe.gov/advanced-analysis/jaws` and cd into jaws/examples/subworkflow-helloworld.

Then run

.. code-block:: bash

   jaws run submit main.wdl inputs.json $(pwd)/out nersc
