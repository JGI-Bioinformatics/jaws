#!/bin/bash -l
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -N 1
#SBATCH -C skylake
#SBATCH -A fungalp
#SBATCH -q jgi_shared
#SBATCH -t 00:10:00

sleep 1000
