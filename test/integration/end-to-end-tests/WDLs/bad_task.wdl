workflow fq_count {
    File fastq_file
    call count_seqs { input: infile = fastq_file }
    output {
        File outfile = count_seqs.outfile
    }
}

task count_seqs {
    File infile
    command <<<
        set -euo pipefail
        bad_cmd_name -l ${infile} | perl -ne 'if (/^\s*(\d+)/ and !($1%4)) {print $1/4, " sequences\n"} else {print STDERR "Invalid Fastq file\n"}' > num_seqs.txt
    >>>
    output {
        File outfile = "num_seqs.txt"
    }
    runtime {
        memory: "5G"
        time: "00:10:00"
        cpu: 1
    }
}

