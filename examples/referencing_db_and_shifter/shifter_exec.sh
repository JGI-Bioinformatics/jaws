#!/bin/bash

set -ex
shifter --image=$1 -V /global/dna/shared/rqc/ref_databases/ncbi/CURRENT:/refdata $2 $3

#${job_shell} ${script}
#shifter --image= exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata docker://$3 $4 $5

