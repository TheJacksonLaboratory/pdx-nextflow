process FASTQ_SORT {

  tag "$sampleID"

  cpus 1
  memory { 50.GB * task.attempt }
  time { 2.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/fastq-tools.sif'

  input:
  tuple val(sampleID), file(trimmed_hsa)

  output:
  tuple val(sampleID), file("*sorted_human*{1,2}.fastq"), emit: sorted_fastq

  script:
  log.info "----- Sorting human reads from Xenome step running on: ${sampleID} -----"
  command_two = params.read_type == 'PE' ? "fastq-sort --id ${trimmed_hsa[1]} > ${sampleID}_sorted_human_2.fastq" : ''

  """
  fastq-sort --id ${trimmed_hsa[0]} > ${sampleID}_sorted_human_1.fastq
  ${command_two}
  """
}
