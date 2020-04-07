#!/bin/bash
shifter --image=$1 -V /global/dna/shared/rqc/ref_databases:/refdata $2 $3
