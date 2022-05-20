====================
WHY Use JAWS
====================

.. role:: bash(code)
  :language: bash

* With a simple command, :bash:`jaws submit`, you can submit a workflow. Running pipelines is made easier by conveniently masking the details of how and where the jobs run.

* You can monitor your jobs with simple commands, like :bash:`jaws log`. 

* You have access to multiple compute sites. For example, to submit to Cori, the command would be :bash:`jaws submit <wdl> <inputs.json> <compute-site>`.
 
* Your workflows are easily scalable. For example, if you have tasks that can be run in parallel, you can select multiple nodes and multiple workers per node.
 
* When developing a pipeline for JAWS, the compute resources can be specified to match the requirements of each task. Therefore, any user of the pipeline can assume that the pipeline will run efficiently.


Here are some additional benefits to using JAWS
-----------------------------------------------

.. figure:: /Figures/why_jaws.svg
   :scale: 50%



Some technologies used
----------------------

.. figure:: /Figures/logos.png
   :scale: 50%


