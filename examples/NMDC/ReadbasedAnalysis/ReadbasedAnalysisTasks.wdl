task profilerGottcha2 {
    Array[File] READS
    String OUTPATH
    String PREFIX
    String? RELABD_COL = "ROLLUP_DOC"
    String DOCKER
    Int? CPU = 4

    command <<<
        mkdir -p ${OUTPATH}

        gottcha2.py -r ${RELABD_COL} \
                    -i ${sep=' ' READS} \
                    -t ${CPU} \
                    -o ${OUTPATH} \
                    -p ${PREFIX} \
                    --database /refdata/nmdc/gottcha2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna
        
        grep "^species" ${OUTPATH}/${PREFIX}.tsv | ktImportTaxonomy -t 3 -m 9 -o ${OUTPATH}/${PREFIX}.krona.html -
    >>>
    output {
        File orig_out_tsv = "${OUTPATH}/${PREFIX}.full.tsv"
        File orig_rep_tsv = "${OUTPATH}/${PREFIX}.tsv"
        File krona_html = "${OUTPATH}/${PREFIX}.krona.html"
    }
    runtime {
        docker: DOCKER
        cpu: CPU
		poolname: "readbaseanalysis-pool"
		node: 1
		nwpn: 1
		mem: "50G"
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
    String OUTPATH
    String PREFIX
    Int? CPU = 4
    String DOCKER

    command <<<
        mkdir -p ${OUTPATH}

        centrifuge -x /refdata/nmdc/centrifuge/p_compressed \
                   -p ${CPU} \
                   -U ${sep=',' READS} \
                   -S ${OUTPATH}/${PREFIX}.classification.csv \
                   --report-file ${OUTPATH}/${PREFIX}.report.csv
        
        ktImportTaxonomy -m 4 -t 2 -o ${OUTPATH}/${PREFIX}.krona.html ${OUTPATH}/${PREFIX}.report.csv
    >>>
    output {
        File orig_out_tsv = "${OUTPATH}/${PREFIX}.classification.csv"
        File orig_rep_tsv = "${OUTPATH}/${PREFIX}.report.csv"
        File krona_html = "${OUTPATH}/${PREFIX}.krona.html"
    }
    runtime {
        docker: DOCKER
        cpu: CPU
		poolname: "readbaseanalysis-pool"
		node: 1
		nwpn: 1
		mem: "50G"
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
    String OUTPATH
    String PREFIX
    Boolean? PAIRED = false
    Int? CPU = 4
    String DOCKER

    command <<<
        mkdir -p ${OUTPATH}
        
        kraken2 ${true="--paired" false='' PAIRED} \
                --threads ${CPU} \
                --db /refdata/nmdc/kraken2/ \
                --output ${OUTPATH}/${PREFIX}.classification.csv \
                --report ${OUTPATH}/${PREFIX}.report.csv \
                ${sep=' ' READS}

        ktImportTaxonomy -m 3 -t 5 -o ${OUTPATH}/${PREFIX}.krona.html ${OUTPATH}/${PREFIX}.report.csv
    >>>
    output {
        File orig_out_tsv = "${OUTPATH}/${PREFIX}.classification.csv"
        File orig_rep_tsv = "${OUTPATH}/${PREFIX}.report.csv"
        File krona_html = "${OUTPATH}/${PREFIX}.krona.html"
    }
    runtime {
        docker: DOCKER
        cpu: CPU
		poolname: "readbaseanalysis-pool"
		node: 1
		nwpn: 1
		mem: "50G"
		time: "04:00:00"
		shared: 1
    }
    meta {
        author: "Po-E Li, B10, LANL"
        email: "po-e@lanl.gov"
    }
}
