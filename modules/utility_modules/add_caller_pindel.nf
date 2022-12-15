process ADD_CALLER_PINDEL {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container "/pdx/pdx_resource_service/elion/containers/python_2.7.3.sif"

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'pindel' }", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("${sampleID}_microIndels.DPfiltered.vcf"), emit: vcf

  script:
  
  """
  bash ${projectDir}/bin/exome/caller_add_pindel.sh ${vcf} ${sampleID}_microIndels.DPfiltered.vcf
  """
}
