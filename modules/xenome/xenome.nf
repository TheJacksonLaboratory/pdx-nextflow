process XENOME_CLASSIFY {

  tag "$sampleID"

  cpus 8
  memory { 50.GB * task.attempt }
  time { 2.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container '/pdx/pdx_resource_service/elion/containers/xenome_1.0.1.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'xenome' }", pattern: "*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(trimmed)

  output:
  tuple val(sampleID), file("human*.fastq"), emit: xenome_fastq
  tuple val(sampleID), file("*.txt"), emit: xenome_stats

  script:

  read_input = params.read_type == 'PE' ? "-i ${trimmed[0]} -i ${trimmed[1]}" : "-i ${trimmed[0]}"
  pairs = params.read_type == 'PE' ? "--pairs" : ""

  """
  /xenome-1.0.1-r/xenome classify -T 8 -P ${params.ref_prefix} ${pairs} --host-name mouse --graft-name human ${read_input} > ${sampleID}_xenome_stats.txt
  """
}
