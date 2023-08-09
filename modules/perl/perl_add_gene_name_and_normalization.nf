process ADD_GENE_NAME_NORM {
  tag "$sampleID"

  cpus 1
  memory { 5.GB * task.attempt }
  time { 2.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container 'quay.io/jaxpdx/r_perl:latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'rsem' }", pattern: "*Normalized", mode:'copy'

  input:
  tuple val(sampleID), file(genes), file(isoforms)

  output:
  file "*withGeneName"
  file "*.Normalized"
  tuple val(sampleID), file("*genes.results.Normalized"), emit: norm_gene_results

  script:
  
  """
  perl ${projectDir}/bin/rnaseq/GeneName_and_Normalization_without_UCSC.pl \
  -i1 ${genes} \
  -i2 ${isoforms} \
  -a1 ${params.hsa_accession}
  """
}
