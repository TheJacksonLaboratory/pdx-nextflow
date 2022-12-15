process CONCATENATE_READS_PE {

  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '03:00:00'

  container "/pdx/pdx_resource_service/elion/containers/python_2.7.3.sif"

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/concatenated_reads' : 'concatenated_reads' }", pattern: "*fastq.gz", mode:'copy'

  input:
  tuple val(sampleID), file(R1), file(R2)

  output:
  tuple val(sampleID), file("*fastq.gz"), emit: concat_fastq

  script:
  
  """
  cat $R1 > ${sampleID}_R1.fastq.gz
  cat $R2 > ${sampleID}_R2.fastq.gz
  """
}
