================================
Best Practices for Creating WDLs
================================

-----------------------------------
#. Container should have all code, but not necessarily reference data which may be mounted
#. Update version tag whenever the container is changed
#. Reference/use containers by specifying a version, rather than LATEST
#. Conda environments are useful for creating containers; be sure to install by specifying versions explicitly rather than using LATEST
#. Archive containers rather than rely upon a public docker hub, in case need to use an old container to reproduce an analysis
#. Avoid hardcoded paths; include code to fetch reference files/dbs (e.g. ftp)
#. Method should be described, perhaps in a README file in the container or in a separate code repository (e.g. git)
    a. It should have sufficient detail to be comprehensible to a relatively naive user
    b. Any papers, test data, test results should be referenced and/or provided

Workflows/Tasks:
----------------
#. Workflows should be comprised of distinct tasks than can be executed outside the context of the workflow
#. Tasks should produce one/few files, written to current working dir, rather than some user defined directory structure; files can be copied anywhere after the pipeline is done, or as a final task.
    a. Similarly, tasks should take file paths as inputs, not folders
#. Separate tasks which use varying compute resources whenever possible
#. Task requirements should be explicitly document and not assumed
#. Minimize unnecessary assumptions about input to facilitate reuse
#. Grouping tasks into sub-workflows makes them more understandable, reuseable, and more easily maintained
#. Setup sub-workflows/microservices for common, important functions to consolidate code
#. Compute resources should be specified/described for each task (e.g. is it multithreaded?  How many threads recommended?  What are typical RAM requirements?  What is typical runtime?)

Developer community
-------------------
#. Share lessons learned in persistent format (e.g. on wiki, forum; not meeting/presentation)
#. Share and announce your contributions (e.g. git, forum)
#. Favor reuse and contribute to shared tools rather than creating your own or forking

Testing/QA
----------
#. Regression test prior to releasing new container
#. Participate in code reviews with your peers

FAIR data practices
-------------------

.. figure:: /Figures/fair.png
  :scale: 100%


#. Findable:
    a. Metadata and results should be indexed/stored in a data warehouse
#. Accessible:
    a. Metadata and files should be easily retrieved programmatically, e.g. FTP or web service rather than websites. i.e. allow bulk downloads, not just browsing capability
#. Interoperable:
    a. Use common file formats; when homegrown formats are used, they should be generated on the fly as temporary intermediate files
#. Reusability: 
	a. write once, run anywhere
#. Reproducible:
    a. Version workflows and use versioned containers
    b. Release workflows/code in additoin to describing method in paper
    c. Workflow output should contain log indicating:
        #. Workflow and version
        #. Input links to metadata (not just paths); including reference db (metadata must include release version)
        #. Runtime parameters
        #. Auto-generated parameters (e.g. random seeds)
        #. Outputs links to data warehouses (e.g. NCBI IDs, etc.)

