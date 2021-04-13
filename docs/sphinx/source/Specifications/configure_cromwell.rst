###########################################
Configure Backends
###########################################

One thing you can customize in a `Cromwell <https://cromwell.readthedocs.io/en/stable/>`_ configuration file is the backend which, among other things, specifies where your jobs will run. You can, for example, choose to use SLURM, JTM, or your local machine. 

.. note:: **Running in JAWS**

      when running in JAWS, the JTM backend is set for you by default and can't be changed. 

This page describes how you can change backend parameters in the config file for testing your WDL.  When running JAWS you cannot change the configuration settings at all.

*****************************************
How Cromwell uses the Configuration file
*****************************************
Cromwell is installed with a complete default config file. However, for testing, you can overwrite parts of it by specifying another config file that has, for example, docker backend specifications.  By including the -Dconfig.file=<your.config> option to your Cromwell command you can overwrite the default config.  

Example configs

* complete config:  `Cromwell.examples.conf <https://github.com/broadinstitute/Cromwell/blob/develop/Cromwell.example.backends/Cromwell.examples.conf>`_  

* slurm config: `slurm.conf <https://github.com/broadinstitute/Cromwell/blob/develop/Cromwell.example.backends/slurm.conf>`_   

* all official backend examples from Cromwell docs: `all examples <https://github.com/broadinstitute/Cromwell/tree/develop/Cromwell.example.backends>`_

* backend example using: :ref:`singularity_backend`.

* backend example using: :ref:`shifter_backend`.

******************************************
Example of Running a WDL with Slurm Config
******************************************
You can see a working example of how to run a task on slurm.
`working example using slurm <https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples/using_slurm_and_local>`_

Follow the README.md in the above repository to run the WDL (`test.wdl`) and inspect the WDL and `cori.conf` file. 

Essentially, you specify that you want to use slurm in the runtime section of your WDL. If you don't specify anything for the `backend:` then the default will run jobs locally, as set in the config (see below).

.. code-block:: text

	runtime {	
		backend: "SLURM"
	}

And in the config file you would define `SLURM` to be one of the backend `providers`.

.. code-block:: text

	backend {
  		default = "Local"
  		providers {
   			SLURM {
				submit = "sbatch <command>"
			}
	}


