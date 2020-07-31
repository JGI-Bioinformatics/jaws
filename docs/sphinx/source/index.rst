======================
Documentation for JAWS
======================

.. figure:: /Figures/JAWS5.tiny.png
    :scale: 100%

The `Joint Genome Institute <https://jgi.doe.gov>`_  has developed JAWS (JGI's Analysis Workflow Service) as a framework to run
computational workflows. Its purpose is to improve the re-usability and robustness of workflows in an evolving high performance computing `(HPC defined) <https://www.nics.tennessee.edu/computing-resources/what-is-hpc>`_ environment.


.. topic:: Intro

   * JAWS simplifies resource configuration and scheduling
   * Create new pipelines with easy-to-understand `Workflow Definition Language(WDL) <https://software.broadinstitute.org/wdl/>`_

.. toctree::
   :name: extradoc
   :maxdepth: 1
   :caption: Intro

   Why use JAWS <Intro/why_jaws>
   How JAWS Works<Intro/how_jaws>
   FAQ <Intro/faq>
   Contact Us <Intro/contact_us>

.. toctree::
   :name: masterdoc
   :maxdepth: 2
   :caption: Run a Workflow in JAWS

   Quickstart Example <Tutorials/jaws_quickstart>
   JAWS Commands <Tutorials/jaws_usage>
   Defining the Input Data <Tutorials/inputs_json>

.. toctree::
   :name: masterdoc
   :maxdepth: 2
   :caption: Developing your own WDLs

   Conda for WDL Development <Tutorials/wdl_development>
   Docker Images Part I <Tutorials/create_env>
   Docker Images Part II <Tutorials/running_env>
   Write a WDL <Tutorials/building_wdls>
   Registering a Workflow in JAWS <Tutorials/register_wdl_tmp>
   
.. toctree::
   :name: specsdoc
   :maxdepth: 1
   :caption: Configuring Workflows

   Runtime Options <Specifications/configuringJTM_in_wdls>
   Configuring Backends <Specifications/configure_cromwell>

.. toctree::
   :name: specsdoc
   :maxdepth: 1
   :caption: Helpful Resources

   Running Examples <Tutorials/example_wdls>
   Subworkflows <Tutorials/subworkflows>
   Best Practices (WDLs) <Intro/best_practices>
   Current Workflows </Specifications/current_workflows>

.. toctree::
   :name: legaldoc
   :maxdepth: 1
   :caption: Legal

   License <Legal/license>
   Contributors <Legal/contributors>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

