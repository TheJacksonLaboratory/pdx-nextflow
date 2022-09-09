process CLASSIFIER_COVERAGE {
  tag "$sampleID"

  cpus 1
  memory { 15.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/python_2.7.v2.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.bed", mode:'copy'

  input:
  tuple val(sampleID), file(norm_genes)

  output:
  tuple val(sampleID), file("*classification"), emit: classified

  script:
  log.info "----- Classifier and coverage Running on: ${sampleID} -----"

  """
  python ${projectDir}/bin/rnaseq/lymphoma_classifier.py -o ${sampleID}_classification ${norm_genes} ${params.classifier_table} ${sampleID}
  """
}
