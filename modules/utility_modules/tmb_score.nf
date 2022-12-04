process TMB_SCORE {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/bedtools_2.27.1_java_1.8_snpeff_4.3_R.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'tmb' }", mode:'copy'

  input:
    tuple val(sampleID), file(tab), file(count2), file(count3)

  output:
    tuple val(sampleID), file("*HM.tab"), emit: tab
    tuple val(sampleID), file("*score"), emit: score
  
  script:
  if (params.workflow == "ctp")
    """
    bedtools coverage -a ${params.bins_ctpcoverage} -b ${sampleID}_count2 | cut -f 1-5 >> ${sampleID}_HM.tab

    Rscript ${projectDir}/bin/wes/TMB_final_CTP.R ${sampleID}_HM.tab ${sampleID}_TMB.score
    """

  else if (params.workflow == "wes")
    """

    bedtools coverage -a ${params.bins_hexcoverage} -b ${sampleID}_count3 | cut -f 1-5 >> ${sampleID}_HM.tab

    Rscript ${projectDir}/bin/wes/TMB_final_HEX.R ${sampleID}_HM.tab ${sampleID}_TMB.score
    """
}