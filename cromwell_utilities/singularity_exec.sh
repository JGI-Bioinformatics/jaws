#!/bin/bash
export SINGULARITY_CACHEDIR=/global/scratch/jaws/sif_files
export SINGULARITY_PULLFOLDER=/global/scratch/jaws/sif_files
export SINGULARITY_TMPDIR=/global/scratch/jaws/sif_files
export SINGULARITY_LOCALCACHEDIR=/global/scratch/jaws/sif_files
singularity exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata docker://$3 $4 $5
exit $?
