process ADD_CALLER_PINDEL {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'


  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*DPfiltered.vcf"), emit: vcf

  script:
  log.info "----- Computing locus depth for: ${sampleID} -----"

  """
  ${projectDir}/bin/wes/caller_add_pindel.sh ${vcf} ${sampleID}_microIndels.DPfiltered.vcf
  """
}
