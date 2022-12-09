process ADD_CALLER_PINDEL {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'pindel' }", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("${sampleID}_microIndels.DPfiltered.vcf"), emit: vcf

  script:
  
  """
  ${projectDir}/bin/exome/caller_add_pindel.sh ${vcf} ${sampleID}_microIndels.DPfiltered.vcf
  """
}
