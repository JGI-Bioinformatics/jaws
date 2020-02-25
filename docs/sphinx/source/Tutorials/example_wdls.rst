============
Example WDLs
============

*******
Summary
*******

You can run WDLs either in JAWS or using a cromwell.jar installation which is usefull for testing. I have examples for both methods below.

*******************************************
Running WDLs with cromwell.jar installation
*******************************************

These examples will run using either our installation or your installation of cromwell.jar.  These are not meant to be run in JAWS.

* `using SLURM v.s. Local: <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/using_slurm_and_local>`_  

	example of a wdl using slurm and local machines

* `subworkflows: <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/subworkflows_and_conditionals>`_ 

	example of using sub-workflows

* `scattering: <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/scatter_gather_example>`_ 

	example of running jobs in parallel.

* `sharding input files: <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/jaws-sharding>`_ 

	example of splitting a file like fasta, fastq, etc. so you can process them in parallel.

* `using reference databases at cori: <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/referencing_db_and_shifter>`_ 

	example of how to use a nt database for blast on cori.

* `using reference databases at LAB IT (Lawrencium): <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/referencing_db_and_singularity>`_ 

	example of how to use a nt database for blast on singularity.

*****************************
Example WDLs that run on JAWS
*****************************
You need to have the jaws environment activated. See :doc:`Set up the JAWS environment </Tutorials/setting_up_environment>`

Before you run the WDL on JAWS, you need to insert real paths for the inputs data (inputs.align.json) and make sure that jaws has permissions to access those data (world read & execute access for each directory leading up to the data). The example's README will tell you how to run it. 

* `alignment example on JAWS: <https://gitlab.com/jfroula/jaws-example-wdl/tree/master/jaws-alignment-example>`_ 

	example running a pipeline in JAWS. It demonstrates the usage of a docker image, a subworkflow and how sharding of an input file is done.


****************
Additional links
****************
* `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.
* :doc:`Best Practices when Writing WDLs </Intro/best_practices>`
