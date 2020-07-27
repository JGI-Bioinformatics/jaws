#!/bin/bash

if [[ ! -e "initial" ]]; then
	echo there is no file \"initial\"
	echo "create it by command: squeue -u jaws | awk '{print \$1}'|tail -n +2|sort -n > initial"
	exit 1
fi

while :; do
	squeue -u jaws | awk '{print $1,$5}'|tail -n +2|sort -n | tee status | awk '{print $1}' > next
	list=(`diff -23 next initial | grep '<' | tr -d ' <'`)
	
	if [[ $list ]]; then
        for i in ${list[@]}; do 
      		pth=$(scontrol show job $i | grep Command= | awk -F'=' '{print $2}')
      		name=$(grep jtm-worker $pth | awk '{print $13}')
			s=$(grep $i status)
      		echo $s $name; 
		done
	fi
	echo sleeping for 8
    sleep 8
done

