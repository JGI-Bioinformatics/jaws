import "rqcfilter_jgi.wdl" as qc
import "ReadbasedAnalysis.wdl" as tax
import "metagenome_assembly_and_alignment.wdl" as asm_aln

workflow combined_workflow {

  Array[File] reads_files
  String prefix
  Boolean? paired = false
  Boolean nersc = false
  Map[String, Boolean] enabled_tax_tools

  call qc.jgi_rqcfilter {
    input:
      input_files = reads_files
  }

  call tax.ReadbasedAnalysis {
    input:
      enabled_tools = enabled_tax_tools,
      reads = reads_files,
      prefix = prefix,
      paired = paired
  }

  call asm_aln.metagenome_assembly_and_alignment {
    input:
      input_files = reads_files,
      nersc = nersc
  }
}
