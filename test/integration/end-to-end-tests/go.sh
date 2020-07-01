#!/bin/bash

for i in {1..5}; do

	echo $i > tmp.txt
	sleep 5
	if [[ "$i" -gt 4 ]]; then
		break
	fi
done
