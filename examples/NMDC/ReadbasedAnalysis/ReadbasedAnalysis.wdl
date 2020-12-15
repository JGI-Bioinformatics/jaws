import "ReadbasedAnalysisTasks.wdl" as tp

workflow ReadbasedAnalysis {
    Map[String, Boolean] enabled_tools
    Array[File] reads
    Int? cpu = 32
    String prefix
    String outdir
    Boolean? paired = false
    String docker = "microbiomedata/nmdc_taxa_profilers:1.0.0"

    if (enabled_tools["gottcha2"] == true) {
        call tp.profilerGottcha2 {
            input: READS = reads,
                   PREFIX = prefix,
                   CPU = cpu,
                   DOCKER = docker
        }
    }
    if (enabled_tools["kraken2"] == true) {
        call tp.profilerKraken2 {
            input: READS = reads,
                   PAIRED = paired,
                   PREFIX = prefix,
                   CPU = cpu,
                   DOCKER = docker
        }
    }
    if (enabled_tools["centrifuge"] == true) {
        call tp.profilerCentrifuge {
            input: READS = reads,
                   PREFIX = prefix,
                   CPU = cpu,
                   DOCKER = docker
        }
    }
    meta {
        author: "Po-E Li, B10, LANL"
        email: "po-e@lanl.gov"
        version: "1.0.0"
    }
}
