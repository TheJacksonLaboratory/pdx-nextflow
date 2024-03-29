process CLASSIFIER_COVERAGE {
  tag "$sampleID"

  cpus 1
  memory { 15.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container 'quay.io/jaxpdx/python_2.7.v2:latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'classifier' }", pattern: "*classification", mode:'copy'

  input:
  tuple val(sampleID), file(norm_genes)

  output:
  tuple val(sampleID), file("*classification"), emit: classified

  script:
  

  """
  python ${projectDir}/bin/rnaseq/lymphoma_classifier.py -o ${sampleID}_classification ${norm_genes} ${params.classifier_table} ${sampleID}
  """
}
