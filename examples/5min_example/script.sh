#!/bin/bash
READS=$1
REF=$2

if [[ $READS ]]; then
	# align reads to reference contigs
	bbmap.sh in=$READS ref=$REF out=test.sam

	# create a bam file from alignment
	samtools view -b -F0x4 test.sam | samtools sort - > test.sorted.bam
else
	echo "Checking that bbmap.sh and samtools is available"
	which bbmap.sh
	if [[ $? > 0 ]]; then 
		echo "bbmap.sh not found"	
		exit 1
	fi
	which samtools
	if [[ $? > 0 ]]; then 
		echo "bbmap.sh not found"	
		exit 1
	fi
fi
