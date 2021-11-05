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
    backend: "parsl"
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
    backend: "parsl"
  }
}
