#!/bin/bash
export SINGULARITY_CACHEDIR=/global/scratch/jfroula/JAWS/img-omics-wdl/JAWS2.0/sif_files
export SINGULARITY_PULLFOLDER=/global/scratch/jfroula/JAWS/img-omics-wdl/JAWS2.0/sif_files
export SINGULARITY_TMPDIR=/global/scratch/jfroula/JAWS/img-omics-wdl/JAWS2.0/sif_files
export SINGULARITY_LOCALCACHEDIR=/global/scratch/jfroula/JAWS/img-omics-wdl/JAWS2.0/sif_files
singularity exec --bind $1:$2 --bind /global/scratch/jfroula/jaws-data/img:/refdata docker://$3 $4 $5
#singularity exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata $3 $4 $5
