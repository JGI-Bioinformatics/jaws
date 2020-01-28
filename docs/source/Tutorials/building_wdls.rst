=================
How to build WDLs
=================

.. role:: bash(code)
   :language: bash

.. note::
    To start learning WDLs you should follow the official `WDL site <https://software.broadinstitute.org/wdl/documentation/>`_ .
    
Below is a step-by-step example of creating a WDL from bash code and running it in cromwell. This section is actually part of the full tutorial "Add workflows to JAWS: the whole process" 

* start with the official `WDL site <https://software.broadinstitute.org/wdl/documentation/>`_
* link to `real world examples <https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-scripts>`_.


#########################################
Converting Example Bash Code to a WDL
#########################################

If we have a mini bash pipeline
-------------------------------

.. code-block:: bash

   #!/bin/bash

   # align reads to reference contigs
   bbmap.sh in=reads.fq ref=reference.fasta out=test.sam

   # create a bam file from alignment
   samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam


The corresponding WDL would look like
-------------------------------------

.. code-block:: bash

    workflow bbtools {
        File reads
        File ref

        call alignment {
           input: fastq=reads,
                  fasta=ref
        }
        call samtools {
           input: sam=alignment.sam
       }
    }

    task alignment {
        File fastq
        File fasta

        command {
    		shifterimg pull jfroula/bbtools:1.2.1 && \
    		shifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=${fastq} ref=${fasta} out=test.sam
        }
        output {
           File sam = "test.sam"
        }
    }

    task samtools {
    	File sam

        command {
       	    shifter --image=jfroula/bbtools:1.2.1 samtools view -b -F0x4 ${sam} | samtools sort - > test.sorted.bam
        }
        output {
       	    File bam = "test.sorted.bam"
        }
    }


Refer to the official WDL website for deeper description and examples.  I'll just point out quickly what's going on here:

  1) The workflow name is :bash:`bbtools` which is used in the :bash:`inputs.json` file to set :bash:`reads` and :bash:`ref`.

  2) The WDL calls two functions or tasks.  The second task, :bash:`samtools` uses the output from the previous task, :bash:`alignment`.

  3) How to pass the output of one task as input to another:  In this example, each of the two tasks has an output section that defines the name of the output.  The name of the output for the alignment task is "sam" (e.g. :bash:`File sam = \"test.sam\"`). Now the second task :bash:`samtools` can access this output by refering to it as "alignment.sam" (<task><dot><output variable>). See the line :bash:`input: sam=alignment.sam`.

  5) Note that each command, in the "command" stanza, is run in a docker container using shifter.


The input file ("inputs.json") would look like this
---------------------------------------------------

.. code-block:: bash

   {
    "bbtools.reads": "<full_path>/reads.fq",
    "bbtools.ref": "<full_path>/reference.fasta"
   }

Test the WDL
------------
Test by following directions in :doc:`Running Workflows </Tutorials/jaws_usage>`
