workflow bbtools { File reads
    File ref

    call alignment {
       input: fastq=reads,
              fasta=ref
    }
    call samtools {
       input: sam=alignment.sam
   }
}

task alignment {
    File fastq
    File fasta

    command <<<
        shifterimg pull jfroula/bbtools:1.2.1 && \
        shifter --image=jfroula/bbtools:1.2.1 bbmap.sh in=${fastq} ref=${fasta} out=test.sam
    >>>
    output {
       File sam = "test.sam"
    }
}

task samtools {
    File sam

    command {
       shifter --image=jfroula/bbtools:1.2.1 samtools view -b -F0x4 ${sam} | shifter --image=jfroula/bbtools:1.2.1 samtools sort - > test.sorted.bam
    }
    output {
       File bam = "test.sorted.bam"
    }
	#runtime {
	#  docker: "jfroula/bbtools:1.2.1"
	#}
}
