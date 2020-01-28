============
Example WDLs
============

*******
Summary
*******

You can run WDLs either in JAWS or using a cromwell.jar installation directly. I have examples for both methods below.

*******************************************
Running WDLs with cromwell.jar installation
*******************************************

These run using either our installation or your installation of cromwell.jar.  These are not meant to be run in JAWS.

* `slurm_local_example <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/using_slurm_and_local>`_ example of a wdl using slurm and local machines
* `sub-workflows <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/subworkflows_and_conditionals>`_ example of using sub-workflows
* `scattering <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/scatter_gather_example>`_ example of running jobs in parallel.
* `sharding files <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/jaws-sharding>`_ example of splitting a file like fasta, fastq, etc. so you can process them in parallel.

*****************************
Example WDLs that run on JAWS
*****************************

The only required modification to the WDL for it to work on JAWS would be that you include the line **backend: "JTM"** in the runtime{} stanza of the above WDLs.
To start running JAWS, see :doc:`Running Workflows`

.. Warning:: some files have disapeared...not a working example
* `jaws-leo-example <https://gitlab.com/jfroula/jaws-leo-example>`_ is an example running a pipeline in JAWS that has only 3 tasks; it demonstrates how to use the static worker pools ("small", and "medium") as well as a dynamic worker pool(submit to slurm).


****************
Additional links
****************
* `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.
* :doc:`Best Practices when Writing WDLs </Intro/best_practices>`
