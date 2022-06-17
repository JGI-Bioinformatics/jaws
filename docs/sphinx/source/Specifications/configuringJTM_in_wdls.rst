#####################################################
JTM Backend Configuration Options When Creating WDL
#####################################################

.. role:: bash(code)
   :language: bash


Use the following tables to help figure out how to configure your runtime{} section.

How to Allocate Resources in your Runtime Section
-------------------------------------------------
Condor is the back end to Cromwell and is responsible for grabbing the appropriatly sized resource from slurm for each wdl-task.  Condor can determine what resource your task needs from only :bash:`memory` and :bash:`cpu` which is set in the runtime{} section. In fact, :bash:`memory` and :bash:`cpu` have defaults set to "5G" and 2(threads), respectively, so you don't have to include them but it is advised for tranparency sake.

.. note::
	Inside your runtime{} section of the WDL, :bash:`cpu` should be set to threads and not cpus, despite the name, because Condor expects that value to be threads .


****************************
Table of available resources
****************************


+-------------+--------+-------+---------+-----+---------+
|    site     |  pool  | nodes | mem(g)* | hrs | threads |
+=============+========+=======+=========+=====+=========+
| cori(nersc) | medium | 2388  | 118     |  72 |   64    |
+             +--------+-------+---------+-----+---------+
|             | xlarge |  20   | 1450    | 168 |   72    |
+-------------+--------+-------+---------+-----+---------+
| jgi(lab-it) | small  | 316   |  46     |  72 |   64    |
+             +--------+-------+---------+-----+---------+
|             | medium |  72   | 236     |  72 |   64    |
+             +--------+-------+---------+-----+---------+
|             | large  |   8   | 492     |  72 |   64    |
+-------------+--------+-------+---------+-----+---------+
| tahoma(emsl)| medium | 160   | 364     |  48 |   72    |
+             +--------+-------+---------+-----+---------+
|             | xlarge |  24   | 1480    |  48 |   72    |
+-------------+--------+-------+---------+-----+---------+
| AWS         |   (see comment below)                    |
+-------------+--------+-------+---------+-----+---------+

AWS is a valid site for JAWS. However, since it uses it's own scheduling system, simply specify the memory and cpu requirements for each task in the :bash:`runtime` section.


\* This number is the gigabytes you can actually use, minus the overhead. For example, on cori, a "medium" node is advertized at 128G but since we can only use roughly 10% of that because of overhead, we will reserve 118G in our WDL.


Links to documentation about each cluster
-----------------------------------------
* `Cori cluster <https://www.nersc.gov/systems/cori/>`_    
* `Lawrencium cluster <https://it.lbl.gov/service/scienceit/high-performance-computing/lrc/computing-on-lawrencium/>`_  
* `Tahoma cluster <https://www.emsl.pnnl.gov/MSC/UserGuide/tahoma/tahoma_overview.html>`_  
* `Amazon Instance types <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/instance-types.html#AvailableInstanceTypes>`_  


****************
Runtime Examples
****************

For example, if you had 100 jobs that needed to be scattered and they required 10G each to run, your runtime{} section would look like
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: text
	
	runtime {
		memory: "10G"
		cpu: 2  # this is really one core on machines with 2 threads/core
	}

You could run 10 tasks in parallel on one 128G haswell machine on cori (leaving some memory left over as a buffer). When a process finishes a job, it will accept another job until all 100 are completed.

If your task required 8 threads and "5G" RAM.
+++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: text
	
	runtime {
		memory: "5G"
		cpu: 8   # this is really 4 cores on machines with 2 threads/core
	}


For a cori haswell node, you could run 8 tasks in parallel because (64 threads/8 = 8) and we see there is more than 5G of ram per task since (118G/8 = 14G).

TODO (what about this scenario if run on at jgi):
+++++++++++++++++++++++++++++++++++++++++++++++++
If your task required 20G and 2 threads: 

.. code-block:: text
	
	runtime {
		memory: "20G"
		cpu: 2   
	}


You can allocate more than one node by asking for two times the total number of threads.  
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This will give you two cori haswell nodes, but your memory still can't exceed 118G unless you are running an MPI job.

.. code-block:: text
	
	runtime {
		memory: "5G"
		cpu: 128   
	}


