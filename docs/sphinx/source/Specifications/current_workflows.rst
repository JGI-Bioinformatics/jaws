==================
Official JGI WDLs
==================
A few of JGI's workflows have been ported over to WDL. Here you can find a list of them.  

`JGI WDL Pipelines <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories>`_.  

See the tabs 

1. Subgroups and Projects (projects created under this repository).  
2. Shared Projects (projects created under some other repository).  


=======================================
JAWS Maintained Sub-workflows and Tasks
=======================================
These are sub-workflows and tasks that are maintained by the JAWS team that can be used when developing your own WDLs.  

`JGI Subworkflows and Tasks <https://code.jgi.doe.gov/official-jgi-workflows/jgi-wdl-tasks>`_

There are two ways to import a subworkflow into your main.wdl

1. You can import the URL from your main WDL like

.. code-block:: text

	import "https://raw.githubusercontent.com/HumanCellAtlas/skylab/optimus_v2.0.0/library/tasks/FastqToUBam.wdl" as FastqToUBam


2. or using a file path,

.. code-block:: text
	
	import ./FastqToUBam.wdl as FastqToUBam
