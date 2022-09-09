process COVCALC_GATK {
  tag "$sampleID"

  cpus 1
  memory { 15.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/python_2.7.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.bed", mode:'copy'

  input:
  tuple val(sampleID), file(txt)
  val(filename)

  output:
  tuple val(sampleID), file("*.bed"), emit: bed

  script:
  log.info "----- GATK COVCALC Running on: ${sampleID} -----"

  """
  python ${projectDir}/bin/rnaseq/coveragecalculator.py ${txt} ${sampleID}_${filename}_avg_median_coverage.bed
  """
}
