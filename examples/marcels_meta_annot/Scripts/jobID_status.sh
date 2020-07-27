#!/bin/bash
while read i; do idp=$(echo $i|awk '{print $1}');pth=$(scontrol show job $idp | grep Command= | awk -F'=' '{print $2}');name=$(grep jtm-worker $pth | awk '{print $13}'); echo $i $name; done < <(squeue -u jaws|tail -n +2 | awk '{print $1,$5}') 
