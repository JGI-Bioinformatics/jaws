#!/bin/sh

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
##$ -j y -o worker.$TASK_ID.log
#$ -j y -o worker.log
#$ -l normal.c
##$ -l exclusive.c
#$ -N tfmq_blast
##$ -l h_rt=00:30:00
##$ -l ram.c=1G
#$ -pe pe_slots 32

module load python
module unload tfmq
module load tfmq
#module unload tfmq/prod
#module load tfmq/dev
#module load taskfarmermq
#module load taskfarmermq/2.1

#for i in {1..8}
#for i in {1..4}
for i in {1..4}
do
   echo "start worker $i"
   tfmq-worker -q blastmq -t 0 -z &
done
wait
#tfmq-worker -b 5 -t 60 -q test10 

##
## To start the client,
## module load tfmq
## tfmq-client -i taskEnv.lst -w 0 -r 1 -q tfmqtaskqueue
##

