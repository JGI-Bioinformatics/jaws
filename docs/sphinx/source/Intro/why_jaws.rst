#############
WHY Use JAWS
#############

.. role:: bash(code)
  :language: bash

| * With a simple command :bash:`jaws run submit` you can submit a workflow. Running pipelines is made easier by conveniently masking the details of how and where the jobs run.
| 
| * You can monitor your jobs with simple commands, like :bash:`jaws run log`. 
| 
| * You have access to multiple compute sites. For example, to submit to cori, the command would be :bash:`jaws run submit <wdl> <inputs.json> <outdir> cori`  
| 
| * Your workflows are easily scalable. If you have tasks that can be run in parallal, you can select multiple nodes and multiple workers per node. 
| 
| * When developing a pipeline for JAWS, the compute resources can be made to match the requirements of each task. Therefore, any user of the pipeline can assume that the pipeline will run efficiently.

|

Here are some additional benifits to using JAWS
-----------------------------------------------

.. figure:: /Figures/why_jaws.svg
   :scale: 50%

|

Some technologies used
----------------------

.. figure:: /Figures/logos.png
   :scale: 50%

