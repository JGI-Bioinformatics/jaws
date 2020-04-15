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

#### Using sub-workflows
* `subworkflows: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/subworkflows_and_conditionals>`_ 


#### How to run jobs in parallel.
* `scattering: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/scatter_gather_example>`_ 


#### How to split a file (like fasta, fastq, txt, etc.) so you can process each shard in parallel.
* `sharding input files: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/jaws-sharding>`_ 


#### How to use a database (i.e. for blast)
* `using reference databases at cori: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/referencing_db_and_shifter>`_ 


*****************************
Example WDLs that run on JAWS
*****************************
You need to have the jaws environment activated. See :doc:`Set up the JAWS environment </Tutorials/setting_up_environment>`

Before you run the WDL on JAWS, you need to insert real paths for the inputs data (inputs.align.json) and make sure that jaws has permissions to access those data (world read & execute access for each directory leading up to the data). The example's README will tell you how to run it. 

* `alignment example on JAWS: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/jaws-alignment-example>`_ 

	example running a pipeline in JAWS. It demonstrates the usage of a docker image, a subworkflow and how sharding of an input file is done.


****************
Additional links
****************
* `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.
* :doc:`Best Practices when Writing WDLs </Intro/best_practices>`
