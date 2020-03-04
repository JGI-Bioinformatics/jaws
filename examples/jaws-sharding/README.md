# Scripts for Sharding
This repo contains script to perform sharding for WDLs using JAWS. Sharding is the act of splitting up fastq, fasta, etc, files so they can be used in parallel jobs.

# What You Can Shard
* fasta files
* fastq files

# Docker Image
see [jaws-sharding](https://cloud.docker.com/repository/docker/jfroula/jaws-sharding)

Only one docker image required for sharding
```
* jfroula/jaws-sharding:1.0.6
```

Additional image required for this alignment example:
```
* jfroula/aligner-bbmap:1.0.0
```

# Implementing Sharding in Your WDL 
Please see the `shard_test.wdl` 

# Run the Example
```
java -jar /global/dna/projectdirs/DSI/workflows/cromwell/java/cromwell.jar run shard_test.wdl -i shard.json

```


