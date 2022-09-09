process QUALITY_STATISTICS {

  tag "$sampleID"

  cpus 1
  memory { 30.GB * task.attempt }
  time { 4.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/python_2.7.3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'quality_stats' }", pattern: "*.gz_stat", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.gz_stat"), emit: quality_stats
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:
  log.info "----- Quality Stats Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    mode_HQ="-S -M"
    inputfq="${fq_reads[0]}"
  }
  if (params.read_type == "PE"){
    mode_HQ="-M"
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
  }

  """
  python ${projectDir}/bin/shared/filter_trim.py $mode_HQ ${params.min_pct_hq_reads}  $inputfq
  """
}
