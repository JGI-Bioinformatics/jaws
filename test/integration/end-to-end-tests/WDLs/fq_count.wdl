workflow fq_count {

    File fastq_file
    call count_seqs { input: infile = fastq_file }
}

task count_seqs {
    File infile

    parameter_meta {
      infile: "A description of what the fastq should be and what used for."
      mo: "A description of what the fastq should be and what used for."
    } 

    command <<<
        wc -l ${infile} | perl -ne 'if (/^\s*(\d+)/ and !($1%4)) {print $1/4, " sequences\n"} else {print STDERR "Invalid Fastq file\n"}' > num_seqs.txt
        echo image: $SHIFTER_IMAGEREQUEST
    >>>

    output {
        File outfile = "num_seqs.txt"
    }

    runtime {
        poolname: "test_small"
		docker: "ubuntu@sha256:b5a61709a9a44284d88fb12e5c48db0409cfad5b69d4ff8224077c57302df9cf"
		cpu: 1
        node: 1
        nwpn: 1
        memory: "5G"
        time: "00:10:00"
        shared: 0
    }
}

