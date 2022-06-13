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
	The values for :bash:`cpu` are threads and not cores, as you might think, because Condor operates in thread units.


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


****************************
Table of available resources
****************************

Site: CORI
++++++++++
NERSC - JAWS runs on shared nodes (haswell, skylake) and dedicated nodes (genepool) which are part of the `cori cluster <https://www.nersc.gov/systems/cori/>`_.

+----------+-----+------------+----------------+-----+-------+--------------+
|constraint|nodes|  ram (G)   | qos            |cores|threads|max time (hrs)|
+==========+=====+============+================+=====+=======+==============+
| haswell  |2,388|  128 (118*)|genepool_special| 32  |   64  |  72          |
+----------+-----+------------+----------------+-----+-------+--------------+
| haswell  |2,388|  128 (118*)|genepool        | 32  |   64  |  72          |
+----------+-----+------------+----------------+-----+-------+--------------+
|     knl  |9,489|   96 (87*) | regular        | 68  |  272  |  48          |
+----------+-----+------------+----------------+-----+-------+--------------+
| skylake  |  20 |1500 (1450*)| jgi_exvivo     | 36  |   36  | 168          |
+          +     +------------+----------------+-----+-------+--------------+
|          |     |  768 (740*)| jgi_shared     | 18  |   18  | 168          |
+----------+-----+------------+----------------+-----+-------+--------------+

 \* the actual number of gigabytes you can use to reserve memory space for system processes.

|

Site: JGI
+++++++++
LAB IT - JAWS runs on shared nodes (lr3) and dedicated nodes (jgi) which is part of the `lawrencium cluster <https://it.lbl.gov/service/scienceit/high-performance-computing/lrc/computing-on-lawrencium/>`_

+---------+------------------+-----+----------+-----+-------+--------------+
|partition|    constraint    |nodes| ram (G)  |cores|threads|max time (hrs)|
+=========+==================+=====+==========+=====+=======+==============+
|     lr3 |                  | 316 |  64 (45*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+
|     lr3 | lr3_c32,jgi_m256 |  32 |256 (236*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+
|     lr3 | lr3_c32,jgi_m512 |   8 |512 (492*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+
|     jgi |                  |  40 |256 (236*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+

\* the actual number of gigabytes you can use to reserve memory space for system processes.

|

Site: TAHOMA
++++++++++++
Pacific Northwest National Labs -  JAWS runs on dedicated (mscjgi) nodes which is on the `tahoma cluster <https://www.emsl.pnnl.gov/MSC/UserGuide/tahoma/tahoma_overview.html>`_

+----------+------------------+-----+------------+-----+-------+--------------+
|partition |    constraint    |nodes| ram (G)    |cores|threads|max time (hrs)|
+==========+==================+=====+============+=====+=======+==============+
|          |                  | 160 |  384 (364*)|  36 |  36   |      72      |
+----------+------------------+-----+------------+-----+-------+--------------+
| analysis |                  |  24 |1500 (1480*)|  36 |  36   |      72      |
+----------+------------------+-----+------------+-----+-------+--------------+

\* the actual number of gigabytes you can use to reserve memory space for system processes.


.. raw:: html

  <details>
  <summary><a>Example of requesting high-mem nodes from Tahoma</a></summary>

  Using 1500G memory machines
  <br>
  <code>
	<pre>
    runtime {
      partition: "analysis"
      time: "00:30:00"
      memory: "1480G"
      poolname: "highmem_test"
      node: 1
      nwpn: 1
    }
  </pre>
  </code>
  </details>

|

Amazon Web Services (AWS)
+++++++++++++++++++++++++
JAWS runs on AWS
cluster: `Instance types <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/instance-types.html#AvailableInstanceTypes>`_

.. raw:: html

  <details>
  <summary><a>Example of requesting resources for AWS</a></summary>

  <br>
  <code>
	<pre>
    runtime {
      memory: "118G"
      cpu: 16
    }
  </pre>
  </code>
  </details>

|

.. _requesting-workers:

