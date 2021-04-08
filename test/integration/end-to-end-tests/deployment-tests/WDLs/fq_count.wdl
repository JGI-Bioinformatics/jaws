version 1.0

workflow fq_count {
	input {
	  File fastq_file
	}

    call count_seqs { input: infile = fastq_file }

    output {
        File outfile = count_seqs.outfile
    }
}

task count_seqs {
	input {
	  File infile
	}

	meta {
	  volatile: true
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
        time: "00:30:00"
        shared: 0
    }
}

