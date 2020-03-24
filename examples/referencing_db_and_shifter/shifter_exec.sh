#!/bin/bash

set -ex
shifter --image=$1 -V /global/dna/shared/rqc/ref_databases/ncbi/CURRENT:/refdata $2 $3

<<<<<<< HEAD
=======
#${job_shell} ${script}
#shifter --image= exec --bind $1:$2 --bind /global/scratch/jaws/ref_data:/refdata docker://$3 $4 $5

>>>>>>> 8fdd6c85f29e50f51ba78169b488b7c9a681c11d
