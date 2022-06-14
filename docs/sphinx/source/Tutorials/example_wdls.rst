================
Working Examples
================
All these working examples demonstrate prefered ways to solve some of the common problems encountered when developing WDLs. These examples are not really meaningful for non-WDL developers.

.. warning::
    Remember to have the jaws environment activated. See :doc:`Set up JAWS Environment </Tutorials/jaws_quickstart>`


How to Run an Example
---------------------
To run the following examples, click on an example link and clone the repo. For the subworkflow example, you would do:

.. code-block:: text

    git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git

    # If you clicked on the subworkflow link below, you can see the directory name is subworkflows_and_conditionals 
    cd subworkflows_and_conditionals

And then follow the command in the README.md of that example.

.. code-block:: text
    
    jaws submit main.wdl inputs.json cori


Example List
------------

Subworkflows promote the reuse of workflows and help modularize workflows for easier maintainability. 

    * `subworkflows: <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples/tree/master/subworkflows_and_conditionals>`_   

    * `alignment example: <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples/tree/master/jaws-alignment-example>`_ 
    

Scattering is the prefered way to run jobs in parallel. "Scatter-gathering" represents the parallization of jobs and the subsequent combining of results into one array for downstream usage.

    * `scattering: <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples/tree/master/scatter_gather_example>`_ 
    

"Sharding" is our term for splitting a file (like fasta, fastq, txt, etc.) into pieces that can then be processed in parallel.  

    * `sharding input files: <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples/tree/master/jaws-sharding>`_ 
    
Here is an example of how to use a database (e.g. for blast) within your WDL. If you require reference files for your pipeline, you need to have a JAWS administrator copy them to a directory dedicated to reference files. This path is available to use within the WDL if your commands are running outside a docker container. If your commands are running within a docker container, you would use "/refdata"; this path is consistent between all compute sites making your WDL portable.

    * `using reference databases at cori: <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples/tree/master/referencing_db_and_shifter>`_ 

****************
Additional links
****************
Real world general examples 

    * `biowdl <https://github.com/biowdl>`_

    * `dockstore.org <https://dockstore.org/search?_type=workflow&descriptorType=wdl&descriptorType=WDL&searchMode=files>`_

