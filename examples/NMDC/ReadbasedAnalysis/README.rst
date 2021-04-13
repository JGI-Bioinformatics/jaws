The Read-based Analysis Workflow
================================

Summary
-------

The pipeline takes sequencing files (single- or paired-end) and profiles them using multiple taxonomic classification tools with the Cromwell as the workflow manager.

Workflow Diagram
----------------

.. image:: docs/readbased_analysis_workflow.png
   :align: center
   :scale: 50%

Workflow Dependencies
---------------------

Third party software
~~~~~~~~~~~~~~~~~~~~

- GOTTCHA2: 2.1.6 `(BSD-3-Clause-LANL) <https://github.com/poeli/GOTTCHA2/blob/master/LICENSE>`_
- Kraken2: 2.0.8 `(MIT) <https://github.com/DerrickWood/kraken2/blob/master/LICENSE>`_
- Centrifuge: 1.0.4 `(GPL-3) <https://github.com/DaehwanKimLab/centrifuge/blob/master/LICENSE>`_

Database 
~~~~~~~~

Each profiling tool requires databases stored in sub-directories at `/global/cfs/projectdirs/m3408/aim2/database/`.

- GOTTCHA2 database (gottcha2/): The database `RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna` is built by complete genomes of bacteria, archaea and viruses from RefSeq Release 90.
- Kraken2 database (kraken2/): This is a standard Kraken 2 database, built by NCBI RefSeq genomes.
- Centrifuge database (centrifuge/)

Workflow Availability
---------------------

The workflow is available in GitHub:
https://github.com/microbiomedata/ReadbasedAnalysis

The container is available at Docker Hub (microbiomedata/nmdc_taxa_profilers):
https://hub.docker.com/r/microbiomedata/nmdc_taxa_profilers


Running Workflow in Cromwell
----------------------------

Description of the files:

- `.wdl`: the WDL file for read-based analysis pipeline.
- `.wdl`: the WDL file for tasks of each tool.
- `.json`: the example inputs.json file for the pipeline.
- `.conf`: the conf file for running cromwell.
- `.job`: example sbatch file.

Test datasets
-------------

Zymobiomics mock-community DNA control `(SRR7877884) <https://www.ebi.ac.uk/ena/browser/view/SRR7877884>`_

Inputs
~~~~~~

The input is a json file:
    
- `ReadbasedAnalysis.enabled_tools`: set the value of the tool as `true` to enable different profiling tools
- `ReadbasedAnalysis.db`: specify the path of the database
- `ReadbasedAnalysis.reads`: specify the path of the reads
- `ReadbasedAnalysis.prefix`: specify the prefix of output file names
- `ReadbasedAnalysis.outdir`: specify the path of output directory
- `ReadbasedAnalysis.cpu`: cpu numbers

.. code-block:: JSON

    {
        "ReadbasedAnalysis.enabled_tools": {
            "gottcha2": true,
            "kraken2": true,
            "centrifuge": true
        },
        "ReadbasedAnalysis.db": {
            "gottcha2": "/global/cfs/projectdirs/m3408/aim2/database/gottcha2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna",
            "kraken2": "/global/cfs/projectdirs/m3408/aim2/database/kraken2/",
            "centrifuge": "/global/cfs/projectdirs/m3408/aim2/database/centrifuge/p_compressed"
        },
        "ReadbasedAnalysis.reads": [
            "/path/to/SRR7877884.1.fastq.gz",
            "/path/to/SRR7877884.2.fastq.gz"
        ],
        "ReadbasedAnalysis.paired": true,
        "ReadbasedAnalysis.prefix": "SRR7877884",
        "ReadbasedAnalysis.outdir": "/path/to/ReadbasedAnalysis",
        "ReadbasedAnalysis.cpu": 8
    }


Outputs
~~~~~~~

The workflow creates individual output directories for each tool, including classification results, logs.::

    ReadbasedAnalysis/
    |-- centrifuge
    |   |-- SRR7877884.classification.csv
    |   |-- SRR7877884.krona.html
    |   `-- SRR7877884.report.tsv
    |-- gottcha2
    |   |-- SRR7877884.full.tsv
    |   |-- SRR7877884.krona.html
    |   `-- SRR7877884.tsv
    `-- kraken2
        |-- SRR7877884.classification.csv
        |-- SRR7877884.krona.html
        `-- SRR7877884.report.csv


Requirements for Execution
--------------------------

- Docker or other Container Runtime
- Cromwell or other WDL-capable Workflow Execution Tool
- 50 GB RAM

Version History
---------------

- 0.0.1

Point of contact
----------------

Package maintainer: Po-E Li <po-e@lanl.gov>
