# use bbmap to align a genome or transcriptome reference to a fastq
# v 1.0.0 - 2020-03-02
#   initial version

workflow rqc_alignment {
  File input_fastq
  File genome_reference
  File transcriptome_reference
  
  call run_alignment {
    input: fastq = input_fastq, genome_fasta = genome_reference, transcriptome_fasta = transcriptome_reference
  }

  # output - save it all, output feeds inputs?

}

task run_alignment {
  File fastq
  File genome_fasta
  File transcriptome_fasta

  # outputs
  String align_zip = "align.zip" 

  runtime {
    docker: "bryce911/rqc-pipeline:20200302"
    time: "00:30:00"
    mem: "32GB"
    node: 1
    nwpn: 1
    poolname: "rqc-align"
  }
  
  command {
    cwd=$(pwd)
    /jgi-rqc-pipeline/alignment/alignment.py -f ${fastq} -g ${genome_fasta} -t ${transcriptome_fasta} -nf -o $cwd

    # remove big, unused files
    rm *.subsample.fastq.gz
    rm genome/*shred.fa
    rm transcriptome/*.sam.gz
    zip align.zip * genome/* transcriptome/*
  }
  
  output {
    File align_zip_file = "align.zip"
  }
}


