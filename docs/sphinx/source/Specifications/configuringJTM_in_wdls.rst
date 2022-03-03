###########################################
JTM Backend Configuration Options When Creating WDL
###########################################

.. role:: bash(code)
   :language: bash

.. note:: **worker pools**

   A **worker** is what runs your tasks. You can have many workers on multiple nodes which together is a **worker pool**. The pool size is determined by specifying the number of nodes (node) and number-of-workers-per-node (nwpn).  For example,  :bash:`#workers = node * nwpn`.
   The maximum number for nwpn is 32 on cori since this is the number of threads.  Other machines will have a different number of threads.


****************************
Table of available resources
****************************

Use the following tables to help figure out how many jobs (i.e. workers) you can run per node and how many total nodes you can expect to get.

Cori
----
JAWS runs on a dedicated cluster.

+----------+-----+----------+----------------+-----+-------+--------------+
|constraint|nodes| ram (G)  | qos            |cores|threads|max time (hrs)|
+==========+=====+==========+================+=====+=======+==============+
| haswell  |2,388|128 (118*)|genepool_special| 32  |   64  |  72          |
+----------+-----+----------+----------------+-----+-------+--------------+
| haswell  |2,388|128 (118*)|genepool**      | 32  |   64  |  72          |
+----------+-----+----------+----------------+-----+-------+--------------+
|     knl  |9,489| 96 (87*) | regular        | 68  |  272  |  48          |
+----------+-----+----------+----------------+-----+-------+--------------+
| skylake  |  20 |758 (700*)| jgi_exvivo     | 32  |   32  | 168          |
+----------+-----+----------+----------------+-----+-------+--------------+
| skylake  |  20 |250 (230*)| jgi_shared     | 32  |   32  | 168          |
+----------+-----+----------+----------------+-----+-------+--------------+

 \* the actual number of gigabytes you can use to reserve memory space for system processes.

 \** The "genepool" qos is appropriate to use when you have many jobs to submit(>10)
 and not "genepool_special" which is a priority node.

.. raw:: html

  <details>
  <summary><a>See examples of requesting different resources</a></summary>

  Using 250G machines
  <br>
  <code>
	<pre>
    runtime {
      time: "00:30:00"
      memory: "250G"
      poolname: "some-unique-name"
      node: 1
      nwpn: 1
      constraint: "skylake"
      qos: "jgi_shared"
      account: "fungalp"
      cpu: 12
    }
	</pre>
  </code>


  Using 700G machines
  <br>
  <code>
	<pre>
    runtime {
      time: "00:30:00"
      memory: "700G"
      poolname: "some-unique-name"
      node: 1
      nwpn: 1
      constraint: "skylake"
      qos: "jgi_exvivo"
      account: "fungalp"
    }
	</pre>
  </code>

  Using non-priority queue ("genepool")
  <br>
  <code>
  <pre>
    runtime {
      poolname: "some-unique-name"
      node: 1
      nwpn: 1
      memory: "10G"
      time: "00:10:00"
      qos: "regular"
      account: "m342"
    }
	</pre>
  </code>
  </details>

|

JGI
---
JAWS runs on a dedicated clusters LR3 and JGI

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


.. raw:: html

  <details>
  <summary><a>Example of requesting high-mem nodes from JGI</a></summary>

  Using 256G memory machines on lr3. Just by having a memory setting larger than 64G, 
  you will be sent to the jgi partition with 256G nodes.

  <br>
  <code>
    <pre>
    runtime {
      poolname: "highmem_test"
      time: "00:30:00"
      memory: "2360G"
      node: 1
      nwpn: 1
    }
  </pre>
  </code>

  <summary><a>Example of requesting high-mem 512G nodes from JGI</a></summary>

  Here you need to set some more params

  <br>
  <code>
    <pre>
    runtime {
        poolname: "helloworldtest"
        node: 1
        nwpn: 1
        memory: "500G"
        time: "00:10:00"
        account: "lr_jgicloud"
        qos: "condo_jgicloud"
        partition: "lr3"
    }
  </pre>
  </code>

  </details>

|

Pacific Northwest National Labs
-------------------------------
JAWS runs on the Tahoma
cluster: `PNNL <https://www.emsl.pnnl.gov/MSC/UserGuide/compute_resources/cascade_overview.html>`_

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

.. _requesting-workers:

******************
Requesting workers
******************
You request resources in a similar manner as for sbatch jobs. The default options are shown below.  Remember to include quotes for strings.

.. code-block:: text

   runtime {
       time: "00:30:00"         # up to 72hrs
       memory: "5G"             # you get a exclusive machine no matter what this setting is. You have two choices: ["115G"|"500G"]
       poolname: "small"        # your choice.
       node: 1                  # number of nodes in the pool. You only need to set this higher when you are scattering a job.
       nwpn: 1                  # number of workers per node (max is number of threads).  This depends on the job's memory & thread requirements.
       cpu: 1                   # this is not used by JTM if run on cori. You can ignore this parameter until we add other "cluster" options.
       constraint: "haswell"    # [haswell|knl|skylake]. Don't use constraint at all if you want to use the default haswell nodes.
                                # Warning: using "knl" will limit your pool to the debug queue which is 30min. limit (until further notice).
                                # If you want to use high-mem node, set it as "skylake".
   }

If you wanted to use all defaults, you could get away with just specifying poolname.

.. code-block:: text

   runtime {
        poolname: "my_pool_name"
   }


How to estimate the number of workers you will need
---------------------------------------------------------------
**workers = node x nwpn**

You will only need more than one worker if you are running jobs in parallel (e.g. using the scattering function in your WDL).
Lets say you are scattering 100 jobs, and you decide 10 workers will give you the desired speedup (roughly 10x), how would you configure the "runtime{}" section to get 10 workers?
The answer depends on how much memory and threads each job will take (e.g. jobs may have variable memory usage so take the highest value seen in your testing). This assumes you did some profiling of your code (even if it was using "memtime" to get max memory estimates for a job).

The decision process should go something like this:

  1. Decide if you want a regular machine (128G) or a large memory machine (512G). Remember that there is an overhead of roughly 13G that you need to subtract from the total memory, so you'd use memory: "115G" or memory: "500G".
  2. If your job maximum memory usage was 50G, and you are using a regular 115G machine then you can run 2 jobs per node. To get 10 workers, you would request :bash:`node: 5` and :bash:`nwpn: 2`.
  3. Alternatively, if your job max memory usage is 2G and it only uses 1 thread, then set :bash:`node: 1` and :bash:`nwpn: 56` (equals 112G total ram). Remember that nwpn: 64 is the maximum.


for example:
**scattering high memory jobs**

.. code-block:: text

   runtime {
     poolname: "my_pool_name"
     time: "2:00:00"
     memory: "115G"
     node: 5
     nwpn: 2
   }


How many threads do I get per worker
------------------------------------
The answer is "It depends on how many workers you ask for".  Consider the following:
Assuming we have a node with 64 threads. If you wanted to run `blastn -num_threads 4` in parallel, and if memory was not a bottleneck, you could run up to 16 blast tasks (64/4=16) on one node. This would equate to 16 workers per node.

.. code-block:: text

   runtime {
     node: 1
     nwpn: 16
   }



.. note::
   If you re-use a worker pool (e.g. same poolname), make sure to include all the runtime parameters you used in the initial runtime, for all the runtimes.  Let's say you were to define a pool with various non-default parameters, and then used the same poolname in another task *without* specifying all the initial parameters.  If the pool were to timeout or crash for some reason, the second task would be trying to use a pool that doesn't exist anymore and hang.  So by copying all the same runtime parameters for each task, even if it were to timeout, a new pool would be created and the job will run.


*********************************
Example Cases and Best-practices
*********************************

If you want to scatter a task use a pool of >1 workers. For instance, If you have a hundred scatter jobs, having 10 workers will give you a 10x speedup. You can configure how many workers (jobs) you want on a node; this depends on the memory requirements per job. Assuming here that each job takes max of 20G, you could run a max of 5 jobs per node.

.. code-block:: text

   runtime {
       cluster: "cori"
       time: "1:00:00"
       memory: "115G"
       poolname: "my_pool_name"
       node: 2
       nwpn: 5
   }

To re-use a worker pool, copy all the params, not just the name.  In this example, the first task takes 20 minutes and the second task takes 40 minutes so the total needs to be at least 1hr.

.. code-block:: text

   task trim {
      runtime {
        cluster: "cori"
        time: "1:00:00"
        memory: "115G"
        poolname: "my_pool_name"
        node: 1
        nwpn: 10
      }
   }
   task assembly {
      runtime {
        cluster: "cori"
        time: "1:00:00"
        memory: "115G"
        poolname: "my_pool_name"
        node: 1
        nwpn: 10
      }

   }

