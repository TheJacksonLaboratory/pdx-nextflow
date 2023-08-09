process ADD_CALLER_GATK {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container "quay.io/jaxpdx/python_2.7.3:latest"

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*DPfiltered.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(add_filters)

  output:
  tuple val(sampleID), file("*DPfiltered.vcf"), emit: vcf

  script:

  """
  bash ${projectDir}/bin/exome/caller_add_gatk.sh ${add_filters} ${sampleID}_variants.DPfiltered.vcf
  """
}
