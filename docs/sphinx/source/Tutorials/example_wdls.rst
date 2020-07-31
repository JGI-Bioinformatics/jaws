============
Example WDLs
============

*******
Summary
*******

You can run WDLs either in JAWS or using a cromwell.jar installation which is usefull for testing. I have examples for both methods below.

****************
Working Examples
****************

To run the examples in JAWS, you need to have the jaws environment activated. See :doc:`Set up JAWS Environment </Tutorials/jaws_quickstart>`

Examples

* Using sub-workflows  

    `subworkflows: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/subworkflows_and_conditionals>`_   

    `alignment example: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/jaws-alignment-example>`_ 
    
* How to run jobs in parallel.  

    `scattering: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/scatter_gather_example>`_ 
    
* How to split a file (like fasta, fastq, txt, etc.) so you can process each shard in parallel.  

    `sharding input files: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/jaws-sharding>`_ 
    
* How to use a database (i.e. for blast)  

    `using reference databases at cori: <https://code.jgi.doe.gov/advanced-analysis/jaws/tree/dev/examples/referencing_db_and_shifter>`_ 


****************
Additional links
****************
* Real world examples 

	`biowdl <https://github.com/biowdl>`_

	`dockstore.org <https://dockstore.org/search?_type=workflow&descriptorType=wdl&descriptorType=WDL&searchMode=files>`_

* :doc:`Best Practices when Writing WDLs </Intro/best_practices>`
