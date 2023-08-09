process COVCALC_GATK {
  tag "$sampleID"

  cpus 1
  memory { 15.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container 'quay.io/jaxpdx/python_2.7:latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.bed", mode:'copy'

  input:
  tuple val(sampleID), file(txt)
  val(filename)

  output:
  tuple val(sampleID), file("*.bed"), emit: bed

  script:
  
  """
  python ${projectDir}/bin/rnaseq/coveragecalculator.py ${txt} ${sampleID}_${filename}_interval_avg_median_coverage.bed
  """
}
