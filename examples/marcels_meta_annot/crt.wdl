workflow crt {

  String imgap_input_fasta
  String imgap_project_id
  String container

  call run {
    input:
      input_fasta = imgap_input_fasta,
      project_id = imgap_project_id,
      container=container
  }

  call transform {
    input:
      project_id = imgap_project_id,
      crt_out = run.out,
      container=container
  }

  output {
    File crisprs = transform.crisprs
    File gff = transform.gff
  }
}

task run {

  String jar="java -Xmx1536m -jar /opt/omics/bin/CRT-CLI.jar"
  File   input_fasta
  String project_id
  String container

  command {
    ${jar} ${input_fasta} ${project_id}_crt.out
  }

  runtime {
    time: "1:00:00"
    memory: "86G"
    docker: container
    shared: 1
    memory: "115G"
    poolname: "crt"
    node: 1
    nwpn: 1
  }

  output {
    File out = "${project_id}_crt.out"
  }
}

task transform {

  String jar="java -Xmx1536m -jar /opt/omics/bin/CRT-CLI.jar"
  String transform_bin="/opt/omics/bin/structural_annotation/transform_crt_output.py"
  File   crt_out
  String project_id
  String crt_out_local = basename(crt_out)
  String container

  command {
    set -uo pipefail  # java returns error code 1 even apon success so remove set -e
    mv ${crt_out} ./${crt_out_local}
    tool_and_version=$(${jar} -version | cut -d' ' -f1,6)
    ${transform_bin} ${crt_out_local} "$tool_and_version"
  }

  runtime {
    time: "1:00:00"
    memory: "86G"
    docker: container
    shared: 1
    memory: "115G"
    poolname: "crt"
    node: 1
    nwpn: 1
  }

  output{
    File crisprs = "${project_id}_crt.crisprs"
    File gff = "${project_id}_crt.gff"
  }
}

