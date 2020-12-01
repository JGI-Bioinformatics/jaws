###########################################
JTM Configuration options when creating WDL
###########################################

.. role:: bash(code)
   :language: bash

.. note:: **worker pools**

   A **worker** is what runs your tasks. You can have many workers on multiple nodes which together is a **worker pool**. The pool size is determined by specifying the number of nodes (node) and number-of-workers-per-node (nwpn).  For example,  :bash:`#workers = node * nwpn`.
   The maximum number for nwpn is 32 on cori since this is the number of threads.  Other machines will have a different number of threads.

.. Warning::

    Having multiple workers only makes sense if you are running a task in parallel via "scattering".

****************************
Table of available resources
****************************

Use the following tables to help figure out how many jobs (i.e. workers) you can run per node and how many total nodes you can expect to get.

On Cori, JAWS runs on a dedicated cluster.

+---------+-----+----------+----------------+-----+-------+--------------+
|node type|nodes| ram (G)  | qos            |cores|threads|max time (hrs)|
+=========+=====+==========+================+=====+=======+==============+
| haswell |2,388|128 (118*)|genepool_special| 16  |   32  |  72          |
+---------+-----+----------+----------------+-----+-------+--------------+
|     knl |9,489| 96 (87*) | regular        | 68  |  272  |  48          |
+---------+-----+----------+----------------+-----+-------+--------------+
| skylake |  20 |758 (700*)| jgi_exvivo     | 32  |   32  | 168          |
+---------+-----+----------+----------------+-----+-------+--------------+
| skylake |  20 |250 (230*)| jgi_shared     | 32  |   32  | 168          |
+---------+-----+----------+----------------+-----+-------+--------------+

|

At JGI, JAWS runs on a dedicated clusters LR3 and JGI

+---------+------------------+-----+----------+-----+-------+--------------+
|partition|    constraint    |nodes| ram (G)  |cores|threads|max time (hrs)|
+=========+==================+=====+==========+=====+=======+==============+
|     lr3 |                  | 316 |  64 (45*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+
|     lr3 | lr3_c32,jgi_m256 | 32  |256 (236*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+
|     lr3 | lr3_c32,jgi_m512 | 8   |512 (492*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+
|     jgi |                  | 40  |256 (236*)|  32 |  64   |      72      |
+---------+------------------+-----+----------+-----+-------+--------------+

|


At Pacific Northwest National Labs: `PNNL <https://www.emsl.pnnl.gov/MSC/UserGuide/compute_resources/cascade_overview.html>`_

+-----------+-----+----------+-----+-------+--------------+
|constraint |nodes| ram (G)  |cores|threads|max time (hrs)|
+===========+=====+==========+=====+=======+==============+
| n.a.      | 960 |128 (118*)|  16 |   16  | 168          |
+-----------+-----+----------+-----+-------+--------------+

| * the actual number of gigabytes you should request (remember there is overhead).


.. _requesting-workers:

******************
Requesting workers
******************
You request resources in a similar manner as for sbatch jobs. The default options are shown below.  Remember to include quotes for strings.

.. code-block:: text

   runtime {
       time: "00:30:00"         # up to 72hrs
       mem: "5G"                # you get a exclusive machine no matter what this setting is. You have two choices: ["115G"|"500G"]
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

  1. Decide if you want a regular machine (128G) or a large memory machine (512G). Remember that there is an overhead of roughly 13G that you need to subtract from the total memory, so you'd use mem: "115G" or mem: "500G".
  2. If your job maximum memory usage was 50G, and you are using a regular 115G machine then you can run 2 jobs per node. To get 10 workers, you would request :bash:`node: 5` and :bash:`nwpn: 2`.
  3. Alternatively, if your job max memory usage is 2G and it only uses 1 thread, then set :bash:`node: 1` and :bash:`nwpn: 56` (equals 112G total ram). Remember that nwpn: 64 is the maximum.


for example:
**scattering high memory jobs**

.. code-block:: text

   runtime {
     poolname: "my_pool_name"
     time: "2:00:00"
     mem: "115G"
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
       mem: "115G"
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
        mem: "115G"
        poolname: "my_pool_name"
        node: 1
        nwpn: 10
      }
   }
   task assembly {
      runtime {
        cluster: "cori"
        time: "1:00:00"
        mem: "115G"
        poolname: "my_pool_name"
        node: 1
        nwpn: 10
      }

   }

