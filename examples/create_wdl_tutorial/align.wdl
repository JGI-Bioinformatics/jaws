workflow bbtools { 
	File reads
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

    command {
        bbmap.sh in=${fastq} ref=${fasta} out=test.sam
    }

    output {
       File sam = "test.sam"
    }

    runtime {
        docker: "jfroula/aligner-bbmap:1.1.3"
        poolname: "muysmall"
        shared: 1
        node: 1
        nwpn: 1
        mem: "5G"
        time: "00:10:00"
    }
}

task samtools {
    File sam

    command {
       samtools view -b -F0x4 ${sam} | samtools sort - > test.sorted.bam
    }

    output {
       File bam = "test.sorted.bam"
    }

    runtime {
        docker: "jfroula/aligner-bbmap:1.1.3"
        poolname: "muysmall"
        shared: 1
        node: 1
        nwpn: 1
        mem: "5G"
        time: "00:10:00"
    }
}

