workflow fq_count {

    File fastq_file
    call count_seqs { input: infile = fastq_file }



    output {
        File outfile = count_seqs.outfile
    }
}

task count_seqs {
    File infile

	parameter_meta {
	  infile: "A description of what the fastq should be and what used for."
	  mo: "A description of what the fastq should be and what used for."
	} 

    command <<<
        wc -l ${infile} | perl -ne 'if (/^\s*(\d+)/ and !($1%4)) {print $1/4, " sequences\n"} else {print STDERR "Invalid Fastq file\n"}' > num_seqs.txt
    >>>

    output {
        File outfile = "num_seqs.txt"
    }

    runtime {
        poolname: "test_small"
        node: 1
        nwpn: 1
        mem: "10G"
        time: "00:10:00"
        shared: 0
    }
}

