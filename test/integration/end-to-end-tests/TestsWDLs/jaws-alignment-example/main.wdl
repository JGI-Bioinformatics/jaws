import "alignment.wdl" as align

workflow main_wdl { 
    File fastq
    File reference

    # this task calls the sub-workflow named bbmap_shard_wf which 
    # is the alignment.wdl.  
    # It's output is "merged_bam_file"
    call align.bbmap_shard_wf { 
           input: reads = fastq,
                  reference = reference
    }
    call bam_stats {
           input: merged_bam = bbmap_shard_wf.merged_bam_file
    }
}

task bam_stats {
    String merged_bam

    command {
        reformat.sh in=${merged_bam} out=stdout.fq | \
        bbstats.sh in=stdin.fq out=stats
    }

    output {
        File alignment_stats = "stats"
    }

    runtime {
        docker: "jfroula/aligner-bbmap:2.0.1"
        poolname: "extrasmall"
        shared: 0
        node: 1
        nwpn: 1
        mem: "5G"
        time: "00:20:00"
    }
}

