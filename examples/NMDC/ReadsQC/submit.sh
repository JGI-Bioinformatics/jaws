#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --output=/global/cfs/projectdirs/m3408/aim2/rqc_filter/log/rqc_filter_1393.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 62
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your@mail.com
#SBATCH --constraint=haswell
#SBATCH --account=m3408
#SBATCH --job-name=rqc_1393

java -Dconfig.file=shifter.conf -jar /global/common/software/m3408/cromwell-45.jar run -o options.json -m metadata_out.json -i input.json rqcfilter.wdl
