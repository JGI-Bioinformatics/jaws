import "alignment_sub.wdl" as align

workflow main_wdl_sulsj {
    File fastq
    File reference
    String bbtools_mem

    # this task calls the sub-workflow named bbmap_shard_wf which
    # is the alignment.wdl.
    # It's output is "merged_bam_file"
    call align.bbmap_shard_wf {
           input: reads = fastq,
                  reference = reference,
                  bbtools_mem = bbtools_mem
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
        cpu: 1
        memory: "20G"
    }
}
