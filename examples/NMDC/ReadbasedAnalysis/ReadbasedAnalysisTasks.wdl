task profilerGottcha2 {
    Array[File] READS
    String PREFIX
    String? RELABD_COL = "ROLLUP_DOC"
    String DOCKER
    Int? CPU = 4

    command <<<
        gottcha2.py -r ${RELABD_COL} \
                    -i ${sep=' ' READS} \
                    -t ${CPU} \
                    -o . \
                    -p ${PREFIX} \
                    --database /refdata/nmdc/gottcha2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna
        
        grep "^species" ${PREFIX}.tsv | ktImportTaxonomy -t 3 -m 9 -o ${PREFIX}.krona.html -
    >>>
    output {
        File orig_out_tsv = "${PREFIX}.full.tsv"
        File orig_rep_tsv = "${PREFIX}.tsv"
        File krona_html = "${PREFIX}.krona.html"
    }
    runtime {
        docker: DOCKER
        cpu: CPU
		poolname: "readbaseanalysis-pool"
		node: 1
		nwpn: 1
		memory: "45G"
		time: "04:00:00"
		shared: 1
    }
    meta {
        author: "Po-E Li, B10, LANL"
        email: "po-e@lanl.gov"
    }
}

task profilerCentrifuge {
    Array[File] READS
    String PREFIX
    Int? CPU = 4
    String DOCKER

    command <<<
        centrifuge -x /refdata/nmdc/centrifuge/p_compressed \
                   -p ${CPU} \
                   -U ${sep=',' READS} \
                   -S ${PREFIX}.classification.csv \
                   --report-file ${PREFIX}.report.csv
        
        ktImportTaxonomy -m 4 -t 2 -o ${PREFIX}.krona.html ${PREFIX}.report.csv
    >>>
    output {
        File orig_out_tsv = "${PREFIX}.classification.csv"
        File orig_rep_tsv = "${PREFIX}.report.csv"
        File krona_html = "${PREFIX}.krona.html"
    }
    runtime {
        docker: DOCKER
        cpu: CPU
		poolname: "readbaseanalysis-pool"
		node: 1
		nwpn: 1
		memory: "45G"
		time: "04:00:00"
		shared: 1
    }
    meta {
        author: "Po-E Li, B10, LANL"
        email: "po-e@lanl.gov"
    }
}

task profilerKraken2 {
    Array[File] READS
    String PREFIX
    Boolean? PAIRED = false
    Int? CPU = 4
    String DOCKER

    command <<<
        kraken2 ${true="--paired" false='' PAIRED} \
                --threads ${CPU} \
                --db /refdata/nmdc/kraken2/ \
                --output ${PREFIX}.classification.csv \
                --report ${PREFIX}.report.csv \
                ${sep=' ' READS}

        ktImportTaxonomy -m 3 -t 5 -o ${PREFIX}.krona.html ${PREFIX}.report.csv
    >>>
    output {
        File orig_out_tsv = "${PREFIX}.classification.csv"
        File orig_rep_tsv = "${PREFIX}.report.csv"
        File krona_html = "${PREFIX}.krona.html"
    }
    runtime {
        docker: DOCKER
        cpu: CPU
		poolname: "readbaseanalysis-pool"
		node: 1
		nwpn: 1
		memory: "45G"
		time: "04:00:00"
		shared: 1
    }
    meta {
        author: "Po-E Li, B10, LANL"
        email: "po-e@lanl.gov"
    }
}
