version 1.0

workflow fq_count {
	input {
	  File fastq_file
	  String max_time = "00:10:00"
	  String max_mem = "1G"
	}

    call count_seqs { 
	    input: infile = fastq_file,
		       max_mem = max_mem,
		       max_time = max_time,
	}

    output {
        File outfile = count_seqs.outfile
    }
}

task count_seqs {
    input {
	  File infile
	  String max_time
	  String max_mem
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
        memory: max_mem
        time: max_time
        shared: 0
    }
}

