# IMG Omics Annotation WDL

## Setup

### Requirements

* Cromwell (and associated Java dependencies)
* Input binaries, databases, etc. specified in the `inputs.json` file

### Docker image
This pipeline uses a docker image.  You can see the dockerfile and latest tag(version) here:  
`https://hub.docker.com/repository/docker/bfoster1/img-omics`

### Input

Input options are listed in the `inputs.json` file and include locations of
input databases, locations of necessary executables, and Boolean variables
to toggle run options on and off.

Currently, the inputs have some executables' paths written in using the
absolute system paths on Cori. Additionally, it is assumed that there is a
`bin` directory at the same level as the `wdl` files which contains numerous
scripts for different parts of the annotations. On Cori at NERSC, this
directory is
`/global/dna/projectdirs/microbial/omics/gbp/img-annotation-pipeline/bin` and
can be copied here into this directory to be used in a relative way.

The workflow is written to take in a metagenome or isolate FASTA which is split
into some number of smaller files. The directory structure looks like:

`/path/to/splits/<split number>/file.fna`

The splits should go from 1 to N, and N should be specified in `inputs.json`
as the variable `metagenome_annotation.num_splits`. The workflow assumes that
all FASTA splits have the same file name. So, for example, if we have an
example FASTA `example.fna` split into 3, the directory structure would look
like:

```
/path/splits/1/example.fna
/path/splits/2/example.fna
/path/splits/3/example.fna
```

where "example" should be set as the value of the variable 
`annotation.imgap_project_id` and `/path/splits/` should be set as
the value of the variable `annotation.imgap_input_dir`.

## Instructions

Edit `inputs.json` to pick whichever parts of the annotation pipeline(s) should
be run. Make sure that the locations of the input databases and executables are
correctly set as well.

The annotation workflow is structured:

```
|annotation
|- setup
|- structural annotation
|-- pre-qc
|-- trnascan
|-- rfam
|-- crt
|-- prodigal
|-- genemark
|-- gff_merge
|-- fasta_merge
|-- gff_and_fasta_stats
|-- post-qc
|- functional annotation
|-- ko_ec
|-- smart
|-- cog
|-- tigrfam
|-- superfamily
|-- pfam
|-- cath_funfam
|-- signalp
|-- tmhmm
|-- product_name_assign
```

Run the workflow with the command:

`java -jar <Cromwell> annotation.wdl -i inputs.json`

in this directory.
