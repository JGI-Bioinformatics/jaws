###########################################
JTM Configuration options when creating WDL
###########################################

.. note:: **worker pools**

   A **worker** is what runs your tasks. You can have multiple workers per node which is a **worker pool**. You define how many workers you want by specifying the number of nodes (node) and number-of-workers-per-node (nwpn).  For example,  #workers = node * nwpn.
   The maximum number for nwpn is 64 on cori since this is the number of threads.  Other machines will have a different number of threads.

.. Warning::

	Having multiple workers only makes sense if you are running a task in parallel via "scattering".

******************
Requesting workers
******************
Although there is a wait time because worker pools are allocated by submitting to SLURM, you can re-use pools between tasks and even between workflows.  You request resources in a similar manner as for sbatch jobs. These have a time limit of 72 hours.  The default options are shown below.  Quotes signify strings and should be included.

.. code-block:: bash

   runtime {
       cluster: "cori"          # [cori|lbl]
       time: "12:00:00"         # up to 72hrs
       mem: "115G"              # you get a exclusive machine no matter what this setting is. You have two choices: ["115G"|"500G"]
       poolname: "my_pool_name" # your choice.
       shared: 1                # a setting of 1 will allow other tasks (even from other workflows) to use identical pools if the "poolname" is the same.
                                # a setting of 0 will still allow different tasks within the same WDL to reuse the same "poolname", but prevent any other WDLS from reusing a pool. This guarantees that two identical WDLs running at the same time will be given different worker pools even though the poolname is the same.
       node: 1                  # number of nodes in the pool. You only need to set this higher when you are scattering a job.
       nwpn: 1                  # number of workers per node(up to 32).  This depends on the job's memory & thread requirements.
       cpu: 32                  # this is not used by JTM if run on cori. You can ignore this parameter until we add other "cluster" options.
       constraint: "haswell"    # [haswell|knl|skylake]. Don't use constraint at all if you want to use the default haswell nodes.
                                # Warning: using "knl" will limit your pool to the debug queue which is 30min. limit (until further notice).
                                # If you want to use high-mem node, set it as "skylake".
   }

If you wanted to use all defaults, you could get away with just specifying poolname.  However, this is ill-advised; see "note" under "Re-using Pools".

.. code-block:: bash

   runtime {
        poolname: "my_pool_name"
   }


How to estimate the number of workers you will need
---------------------------------------------------------------
**workers = node x nwpn**

You will only need multiple workers if you are running jobs in parallel (e.g. using the scattering function in your WDL).
Lets say you are scattering 100 jobs, and you decide 10 workers will give you the desired speedup (roughly 10x), how would you configure the "runtime{}" section to get 10 workers?
The answer depends on how much memory and threads each job will take (e.g. jobs may have variable memory usage so take the highest value seen in your testing). This assumes you did some profiling of your code (even if it was using "memtime" to get max memory estimates for a job).

The decision process should go something like this:

  1. decide if you want a regular machine (128G - quicker to schedule) or a xlarge machine (512G - longer to schedule). Remember that there is an overhead of roughly 13G that you need to subtract from the total memory, so you'd use mem: "115G" or mem: "500G".
  2. if your job maximum memory usage was 50G, and you are using a regular 115G machine then you can run 2 jobs per node. You would request node: 5 & nwpn: 2.
  3. if your job max memory usage is 2G and it only uses 1 thread, then set node: 1 & nwpn: 56. Remember that nwpn: 64 is the maximum.


for example:
**scattering high memory jobs**

.. code-block:: bash

   runtime {
     poolname: "my_pool_name"
     cluster: "cori"
     time: "2:00:00"
     mem: "115G"
     node: 5
     nwpn: 2
   }


How many threads do I get per worker
------------------------------------
The answer is "It depends on how many workers you ask for".  Consider the following:
Assuming we have a node with 64 threads. If you wanted to run `blastn -num_threads 4` in parallel, and if memory was not a bottleneck, you could run up to 16 blast tasks (64/4=16) on one node. This would equate to 16 workers per node.

.. code-block:: bash

   runtime {
     node: 1
     nwpn: 16
   }


Re-using Pools
--------------
The advantage of setting "poolname" to some user defined name is that you can re-use the pool for another task that will not have to re-submit to SLURM.  Since the second task is re-using the pool, the time limit must be adequate to run both tasks. As for mem, node and nwpn, remember to set these to the highest number you will encounter in either task. In theory, you could reserve a large machine for a long time and do all tasks on that machine, only having to sbatch once; however, this would circumvent the optimization potential of the workflow engine, which is to pair small tasks with small compute resources.

.. note::
   If you re-use a worker pool (e.g. same poolname), make sure to include all the necessary runtime parameters like cpu, time, etc.  Lets say you define a Dynamic pool as in the above example and then use the same poolname: "my_pool_name" in another task without specifying time, mem, etc.  If the pool were to timeout or crash for some reason, the second task would be trying to use a pool that doesn't exist anymore and hang.  So by copying all the runtime parameters for each task using "my_pool_name", even if it were to timeout, a new pool would be created and the job will run.


*********************************
Example Cases and Best-practices
*********************************

If you want to scatter a task use a pool of >1 workers. For instance, If you have a hundred scatter jobs, having 10 workers will give you a 10x speedup. You can configure how many workers(jobs) you want on a node; this depends on the memory requirements per job. Assuming here that each job takes max of 20G ...

.. code-block:: bash

   runtime {
       cluster: "cori"
       time: "1:00:00"
       mem: "115G"
       poolname: "my_pool_name"
       node: 2
       nwpn: 5
   }

To re-use a worker pool, copy all the params, not just the name.  In this example, the first task takes 20 minutes and the second task takes 40 minutes so the total needs to be at least 1hr.

.. code-block:: bash

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

