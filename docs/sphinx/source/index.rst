======================
Documentation for JAWS
======================

.. role:: chocolate

.. figure:: /Figures/JAWS5.tiny.png
    :scale: 100%

The `Joint Genome Institute <https://jgi.doe.gov>`_ has developed JAWS (JGI Analysis Workflow Service) as a framework to run computational workflows. Its purpose is to improve the re-usability and robustness of workflows in a high performance computing (`HPC <https://www.nics.tennessee.edu/computing-resources/what-is-hpc>`_), cluster, and cloud environments.

Summary of JAWS
---------------

* JAWS should simplify your workflow submissions because it takes care of resource configuration and scheduling

* Before using JAWS, workflows need to be wrapped in the human readable `Workflow Definition Language (WDL) <https://software.broadinstitute.org/wdl/>`_

* JAWS commands can be integrated into larger workflows since everything is `CLI <https://en.wikipedia.org/wiki/Command-line_interface>`_ based.

Getting Started
---------------

* `Create your own WDLs <Tutorials/create_env.html>`_
* `Run an existing WDL in JAWS <Tutorials/jaws_quickstart.html>`_


.. toctree::
   :hidden:
   :name: extradoc
   :maxdepth: 1
   :caption: Intro

   Why use JAWS <Intro/why_jaws>
   What is JAWS <Intro/how_jaws>
   JAWS Video Library<Intro/videos>
   FAQ <Intro/faq>
   Known Issues <Intro/known_issues>
   Getting Help <Intro/contact_us>

.. toctree::
   :hidden:
   :name: getstarted
   :maxdepth: 2
   :caption: Run a Workflow in JAWS

   Quickstart Example <Tutorials/jaws_quickstart>
   JAWS Commands <Tutorials/jaws_usage>
   Defining the Input Data <Tutorials/inputs_json>


.. toctree::
   :hidden:
   :name: wdldevelop
   :maxdepth: 2
   :caption: Developing your own WDLs

   Step 1: Development environment <Tutorials/wdl_development>
   Step 2: Create docker container <Tutorials/create_env>
   Step 3: Write a WDL <Tutorials/building_wdls>

.. toctree::
   :hidden:
   :name: specsdoc
   :maxdepth: 1
   :caption: General Resources

   Best Practices for WDLs <Intro/best_practices>
   Runtime Options <Specifications/configuringJTM_in_wdls>
   Code Snippets  <Tutorials/snippets>
   Registering your Workflow <Tutorials/register_wdl>
   Working Examples <Tutorials/example_wdls>
   Subworkflows <Tutorials/subworkflows>
   Current JGI Workflows </Specifications/current_workflows>
