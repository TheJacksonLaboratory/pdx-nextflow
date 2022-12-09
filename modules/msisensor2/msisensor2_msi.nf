process MSISENSOR2_MSI {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container '/pdx/pdx_resource_service/elion/containers/msisensor2_latest.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'msisensor2' }", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file("*msisensor"), emit: msisensor
  file("${sampleID}_msisensor_dis")
  file("${sampleID}_msisensor_somatic")

  script:
  
  """
  mkdir models

  cp -r ${params.msisensor_model} models

  msisensor2 msi -M models/models_hg38 -t ${bam} -o ${sampleID}_msisensor

  """
}
