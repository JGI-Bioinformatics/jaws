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
   Helpful Hints when Developing WDLs <Intro/helpfull_hints>
   FAQ <Intro/faq>
   Contact Us <Intro/contact_us>

.. toctree::
   :name: masterdoc
   :maxdepth: 2
   :caption: Testing and Development

   Quickstart <Tutorials/wdl_development>
   Write a WDL <Tutorials/building_wdls>
   Set up the JAWS environment <Tutorials/setting_up_environment>

.. toctree::
   :name: masterdoc
   :maxdepth: 2
   :caption: Running in JAWS

   Add workflows to JAWS (the whole process) <Tutorials/whole_wdl_process>
   Running WDLs <Tutorials/jaws_usage>
   Registering a Workflow in JAWS <Tutorials/register_wdl>

.. toctree::
   :name: specsdoc
   :maxdepth: 1
   :caption: Configuring Workflows

   Runtime options <Specifications/configuringJTM_in_wdls>
   Configuring Backends <Specifications/configure_cromwell>
   Using Docker Images </Specifications/running_env>

.. toctree::
   :name: specsdoc
   :maxdepth: 1
   :caption: Helpful Resources

   Examples <Tutorials/example_wdls>
   Subworkflows <Tutorials/subworkflows>
   Best Practices (WDLs) <Intro/best_practices>
   Current Workflows </Specifications/current_workflows>

.. toctree::
   :name: legaldoc
   :maxdepth: 1
   :caption: Legal

   license <Legal/license>
   contributors <Legal/contributors>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

