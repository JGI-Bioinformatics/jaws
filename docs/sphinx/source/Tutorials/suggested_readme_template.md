# JGI Metagenome annotation 

## Summary
This workflow was converted from Torben Nielson's JGI metagenome pipeline.  This workflow performs read filtering and error correction steps before running metaSpades.

## Running Requirements
24hrs with an input of 24G fastq file.

## Input files
expects: fastq, illumina, interleaved, paired-end, >130bp 

## Output files
There should be three output files where `results` contains the metagenome annotations.

```
bbduk.log
spades.log
results
```

## Running Environment
To Recreate the Running Environment

The Docker image and Dockerfile can be found at hub.docker.com.
```
https://cloud.docker.com/u/jfroula/repository/docker/jfroula/jgi_meta_annot_torben
image=jfroula/jgi_meta_annot_torben
```

