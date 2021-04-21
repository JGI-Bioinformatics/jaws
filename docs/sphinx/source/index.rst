======================
Documentation for JAWS
======================

.. role:: chocolate

.. figure:: /Figures/JAWS5.tiny.png
    :scale: 100%

The `Joint Genome Institute <https://jgi.doe.gov>`_ has developed JAWS (JGI Analysis Workflow Service) as a framework to run computational workflows. Its purpose is to improve the re-usability and robustness of workflows in a high performance computing (`HPC <https://www.nics.tennessee.edu/computing-resources/what-is-hpc>`_), cluster, and cloud environments.

.. topic:: In Short...

   * JAWS should simplify your workflow submissions because it takes care of resource configuration and scheduling
   * Before using JAWS, workflows need to be wrapped in an easy-to-understand `Workflow Definition Language (WDL) <https://software.broadinstitute.org/wdl/>`_
   * JAWS presents the user with many commands to manage your jobs. 



Some Definitions
-----------------

+---------------------+------------------------------------------------------------------------------+
|:chocolate:`WDL`     | | The Workflow Description Language is essentially a                         |
|                     | | wrapper around the commands in your pipeline code.                         |
+---------------------+------------------------------------------------------------------------------+
|:chocolate:`Cromwell`| | A workflow engine which takes a WDL and converts it                        |
|                     | | to bash commands that can be run on a "backend".                           |
+---------------------+------------------------------------------------------------------------------+
|:chocolate:`Backend` | | Any compute resource like a private labtop, or a public                    |
|                     | | compute cluster like Cori. In our case we use JTM as a backend.            |
+---------------------+------------------------------------------------------------------------------+
|:chocolate:`JTM`     | | A customized backend that uses different compute sites (i.e. Cori, and JGI)|
|                     | | and is responsible for reserving workers on the cluster.                   |
+---------------------+------------------------------------------------------------------------------+
|:chocolate:`Workers` | | Workers are "processes" that run the Cromwell commands.                    |
|                     | | They can be part of a "worker pool" to process parallel tasks.             |
+---------------------+------------------------------------------------------------------------------+
|:chocolate:`Tasks`   | | Each WDL or run is composed of multiple tasks.                             |
+---------------------+------------------------------------------------------------------------------+


.. toctree::
   :hidden:
   :name: extradoc
   :maxdepth: 1
   :caption: Intro

   Why use JAWS <Intro/why_jaws>
   What is JAWS <Intro/how_jaws>
   JAWS Video Library<Intro/videos>
   FAQ <Intro/faq>
   Contact Us <Intro/contact_us>

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
   :name: specsdoc
   :maxdepth: 1
   :caption: General Resources

   Working Examples <Tutorials/example_wdls>
   Subworkflows <Tutorials/subworkflows>
   Current JGI Workflows </Specifications/current_workflows>
   Using Reference Data </Specifications/refdata>
   Known Issues </Specifications/known_issues>

.. toctree::
   :hidden:
   :name: wdldevelop
   :maxdepth: 2
   :caption: Developing your own WDLs

   Conda for WDL Development <Tutorials/wdl_development>
   Docker Images Part I <Tutorials/create_env>
   Docker Images Part II <Tutorials/running_env>
   Write a WDL <Tutorials/building_wdls>
   Best Practices for WDLs <Intro/best_practices>
   Runtime Options <Specifications/configuringJTM_in_wdls>
   Code Snippets  <Tutorials/snippets>
   Registering your Workflow <Tutorials/register_wdl>

