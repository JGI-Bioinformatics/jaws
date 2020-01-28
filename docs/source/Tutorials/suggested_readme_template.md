# JGI Metagenome annotation 

## Summary
This workflow is based on Torben Nielson's pipeline
This workflow runs some read filtering and error correction steps before running metaSpades.

## The Docker image and Dockerfile can be found here
```
https://cloud.docker.com/u/jfroula/repository/docker/jfroula/jgi_meta_annot_torben
image=jfroula/jgi_meta_annot_torben
```

## Running Requirements
1G for 24hrs

## Input files
expects: fastq, illumina, interleaved, paired-end, >130bp 

## Output files
```
bbduk.log
spades.log
results
```
