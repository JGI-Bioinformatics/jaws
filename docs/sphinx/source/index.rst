======================
Documentation for JAWS
======================

.. figure:: /Figures/JAWS5.tiny.png
    :scale: 100%

The `Joint Genome Institute <https://jgi.doe.gov>`_  has developed JAWS (JGI's Analysis Workflow Service) as a framework to run
computational workflows. Its purpose is to improve the re-usability and robustness of workflows in a high performance computing `(HPC defined) <https://www.nics.tennessee.edu/computing-resources/what-is-hpc>`_ environment.

.. topic:: In Short...

   * JAWS should simplify your workflow submissions because it takes care of resource configuration and scheduling
   * Before using JAWS, workflows need to be wrapped in an easy-to-understand `Workflow Definition Language (WDL) <https://software.broadinstitute.org/wdl/>`_
   * JAWS presents the user with many commands and a dashboard to manage your jobs. 


.. topic:: Some Definitions

	* **WDL** => The Workflow Description Language is essentially a wrapper around the commands in your pipeline code.  
	* **Cromwell** => workflow engine which takes a WDL and converts it to bash commands that can be run on a "backend".  
	* **backend** => any compute resource like your computer, or a public compute cluster. In our case we use JTM as a backend.  
	* **JTM** => customized backend that uses different compute sites (i.e. NERSC, and soon LBNL) and is responsible for reserving workers on the cluster on 
	* **workers** => workers are "processes" that run the `Cromwell <https://cromwell.readthedocs.io/en/stable/>`_ bash commands. They can be part of a "worker pool" to process parallel tasks.  
	* **Tasks** => each WDL or run is composed of multiple tasks. 
	* **Run** => one JAWS submission.  


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
   :name: masterdoc
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

.. toctree::
   :hidden:
   :name: masterdoc
   :maxdepth: 2
   :caption: Developing your own WDLs

   Conda for WDL Development <Tutorials/wdl_development>
   Docker Images Part I <Tutorials/create_env>
   Docker Images Part II <Tutorials/running_env>
   Write a WDL <Tutorials/building_wdls>
   Best Practices for WDLs <Intro/best_practices>
   Runtime Options <Specifications/configuringJTM_in_wdls>
   Registering your Workflow <Tutorials/register_wdl>

.. toctree::
   :hidden:
   :name: legaldoc
   :maxdepth: 1
   :caption: Legal

   License <Legal/license>
   Contributors <Legal/contributors>

