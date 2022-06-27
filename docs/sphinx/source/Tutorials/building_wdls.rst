=================
How to build WDLs
=================

.. role:: bash(code)
   :language: bash

In this tutorial, we will see how a bash wrapper (script.sh) can be made into a WDL.

.. note::
    To really learn WDLs you should follow the official `WDL site <https://software.broadinstitute.org/wdl/documentation/>`_.  However, to get
    and idea of what its all about, continue...


Some useful links

    * start with the official `WDL site <https://software.broadinstitute.org/wdl/documentation/>`_
    * `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.
    * Re-usable subworkflow tasks: `WLD-tasks <https://gitlab.com/jgi-doe/wdl-tasks.git>`_


****************************
Clone the Example Repository
****************************

.. code-block:: text

   git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git
   cd jaws-tutorial-example/5min_example


*************************************
Converting Example Bash Code to a WDL
*************************************

If we have a workflow represented as a script (i.e. script.sh), we can parse it up into WDL tasks.

.. note ::
    Each script you create should execute in and write output to the **current working directory**.

Our script.sh has two steps:

.. code-block:: text

   #!/bin/bash

   # align reads to reference contigs
   bbmap.sh in=reads.fq ref=reference.fasta out=test.sam

   # create a bam file from alignment
   samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam



The corresponding WDL would look like this (call it align_final.wdl).

Notice that I am running all the commands inside a docker container :bash:`jfroula/bbtools:1.2.1`.  The image of which should exist in hub.docker.com and will therefore be pulled localy by cromwell, assuming the backend was configured correctly.

.. code-block:: text

    version 1.0

    workflow bbtools {
        input {
            File reads
            File ref
        }

        call alignment {
           input: fastq=reads,
                  fasta=ref
        }
        call samtools {
           input: sam=alignment.sam
       }
    }

    task alignment {
        input {
            File fastq
            File fasta
        }

        command {
            bbmap.sh in=~{fastq} ref=~{fasta} out=test.sam
        }

        runtime {
            docker: "jfroula/aligner-bbmap:2.0.2"
            time: "00:10:00"
            memory: "5G"
            cpu: 1
        }

        output {
           File sam = "test.sam"
        }
    }

    task samtools {
        input {
            File sam
        }

        command {
           samtools view -b -F0x4 ~{sam} | samtools sort - > test.sorted.bam
        }

        runtime {
            docker: "jfroula/aligner-bbmap:2.0.2"
            time: "00:10:00"
            memory: "5G"
            cpu: 1
        }

        output {
           File bam = "test.sorted.bam"
        }
    }

Refer to the official WDL website for deeper description and examples.  I'll just point out quickly what's going on here:

  1) The workflow name is :bash:`bbtools` which is used in the :bash:`inputs.json` file to set :bash:`reads` and :bash:`ref`.

  2) The WDL calls two functions or tasks.  The second task, :bash:`samtools` uses the output from the previous task, :bash:`alignment`.

  3) How to pass the output of one task as input to another:  In this example, each of the two tasks has an output section that defines the name of the output.  The name of the output for the alignment task is "sam" (e.g. :bash:`File sam = \"test.sam\"`). Now the second task :bash:`samtools` can access this output by refering to it as "alignment.sam" (<task><dot><output variable>). See the line :bash:`input: sam=alignment.sam`.

  5) Note that each command, in the "command" stanza, is run in a docker container. 


The input file ("inputs.json") would look like this
---------------------------------------------------

.. code-block:: text

   {
    "bbtools.reads": "<full_path>/reads.fq",
    "bbtools.ref": "<full_path>/reference.fasta"
   }


An example of running this WDL was described in the last section :ref:`Run with Docker Inside the runtime{} <run with conf>`


Re-use Ready Made Tasks
-----------------------
Check if there are any ready made tasks (as subworkflows) that you can use in your WDL. 

`WDL-tasks <https://gitlab.com/jgi-doe/wdl-tasks.git>`_
