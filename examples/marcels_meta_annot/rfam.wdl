workflow rfam {

  String imgap_input_fasta
  String imgap_project_id
  String imgap_project_type
  Int    additional_threads
  String database_location="/refdata/database"
  String cm="${database_location}"+"Rfam/13.0/Rfam.cm"
  String claninfo_tsv="${database_location}"+"Rfam/13.0/Rfam.claninfo"
  String feature_lookup_tsv="${database_location}"+"Rfam/13.0/Rfam_feature_lookup.tsv"
  String container


  call cmsearch {
    input:
      input_fasta = imgap_input_fasta,
      project_id = imgap_project_id,
      cm = cm,
      threads = additional_threads,
      container=container
  }

  call clan_filter {
    input:
      project_id = imgap_project_id,
      tbl = cmsearch.tbl,
      claninfo_tsv = claninfo_tsv,
      feature_lookup_tsv = feature_lookup_tsv,
      container=container
  }

  call misc_and_regulatory {
    input:
      rfam_gff = clan_filter.rfam_gff,
      project_id = imgap_project_id
  }

  call rrna {
    input:
      rfam_gff = clan_filter.rfam_gff,
      project_id = imgap_project_id
  }

  call ncrna_tmrna {
    input:
      rfam_gff = clan_filter.rfam_gff,
      project_id = imgap_project_id
  }

  output {
    File misc_bind_misc_feature_regulatory_gff = misc_and_regulatory.misc_bind_misc_feature_regulatory_gff
    File rrna_gff = rrna.rrna_gff
    File ncrna_tmrna_gff = ncrna_tmrna.ncrna_tmrna_gff
  }
}

task cmsearch {

  String bin="/opt/omics/bin/cmsearch"

  File   input_fasta
  String project_id
  String cm
  Int    threads
  String container

  command {
    ${bin} --notextw --cut_tc --cpu ${threads} --tblout ${project_id}_rfam.tbl ${cm} ${input_fasta}
  }

  runtime {
    time: "1:00:00"
    memory: "86G"
    docker: container
    shared: 1
    memory: "115G"
    poolname: "rfam"
    node: 1
    nwpn: 1
  }

  output {
    File tbl = "${project_id}_rfam.tbl"
  }
}

task clan_filter {

  String clan_filter_bin="/opt/omics/bin/structural_annotation/rfam_clan_filter.py"
  String project_id
  File   tbl
  String cmsearch_bin="/opt/omics/bin/cmsearch"
  String   claninfo_tsv
  String   feature_lookup_tsv
  String container

  command <<<
    set -euo pipefail
    tool_and_version=$(${cmsearch_bin} -h | grep INFERNAL | perl -pne 's/^.*INFERNAL/INFERNAL/' )
    grep -v '^#' ${tbl} | \
    awk '$17 == "!" {print $1,$3,$4,$6,$7,$8,$9,$10,$11,$15,$16}' | \
    sort -k1,1 -k10,10nr -k11,11n | \
    ${clan_filter_bin} "$tool_and_version" \
    ${claninfo_tsv} ${feature_lookup_tsv} > ${project_id}_rfam.gff
  >>>

  runtime {
    time: "1:00:00"
    memory: "86G"
    docker: container
    shared: 1
    memory: "115G"
    poolname: "rfam"
    node: 1
    nwpn: 1
  }

  output {
    File rfam_gff = "${project_id}_rfam.gff"
  }
}

task misc_and_regulatory {
  
  File   rfam_gff
  String project_id

  command <<<
    awk -F'\t' '$3 == "misc_bind" || $3 == "misc_feature" || $3 == "regulatory" {print $0}' \
    ${rfam_gff} > ${project_id}_rfam_misc_bind_misc_feature_regulatory.gff
  >>>

  runtime {
    time: "1:00:00"
    memory: "86G"
    shared: 1
    memory: "115G"
    poolname: "rfam"
    node: 1
    nwpn: 1
  }

  output {
    File misc_bind_misc_feature_regulatory_gff = "${project_id}_rfam_misc_bind_misc_feature_regulatory.gff"
  }
}

task rrna {

  File   rfam_gff
  String project_id

  command <<<
    awk -F'\t' '$3 == "rRNA" {print $0}' ${rfam_gff} > ${project_id}_rfam_rrna.gff
  >>>

  runtime {
    time: "1:00:00"
    memory: "86G"
    shared: 1
    memory: "115G"
    poolname: "rfam"
    node: 1
    nwpn: 1
  }

  output {
    File rrna_gff = "${project_id}_rfam_rrna.gff"
  }
}

task ncrna_tmrna {

  File   rfam_gff
  String project_id

  command <<<
    awk -F'\t' '$3 == "ncRNA" || $3 == "tmRNA" {print $0}' \
        ${rfam_gff} > ${project_id}_rfam_ncrna_tmrna.gff
  >>>

  runtime {
    time: "1:00:00"
    memory: "86G"
    shared: 1
    memory: "115G"
    poolname: "rfam"
    node: 1
    nwpn: 1
  }

  output {
    File ncrna_tmrna_gff = "${project_id}_rfam_ncrna_tmrna.gff"
  }
}

