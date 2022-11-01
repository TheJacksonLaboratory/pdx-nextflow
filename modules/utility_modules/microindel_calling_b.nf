process MICROINDEL_CALLING_B {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/bedtools_2.27.1_python_2_7_3.sif'


  input:
  tuple val(sampleID), file(all_vcf)

  output:
  tuple val(sampleID), file("*.DPfiltered1.vcf"), emit: vcf

  script:
  log.info "----- Microindel calling, part 2 running on ${sampleID} -----"

  """
  bedtools intersect \
  -header \
  -a ${all_vcf} \
  -b  ${params.targets_gatk} \
  -f 1.0 > ${sampleID}_microIndels.raw.vcf

  python ${projectDir}/bin/wes/filter_for_minimum_depth.py ${sampleID}_microIndels.raw.vcf ${sampleID}_micro_Indels.DPfiltered1.vcf

  touch  ${sampleID}_BP  ${sampleID}_D  ${sampleID}_DSI  ${sampleID}_INT  ${sampleID}_INT_final  ${sampleID}_INV  ${sampleID}_LI  ${sampleID}_RP  ${sampleID}_SI  ${sampleID}_TD


  """
}
