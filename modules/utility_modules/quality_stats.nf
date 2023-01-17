process QUALITY_STATISTICS {

  tag "$sampleID"

  cpus 1
  memory { 30.GB * task.attempt }
  time { 4.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container '/pdx/pdx_resource_service/elion/containers/python_2.7.3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'quality_stats' }", pattern: "*_fastqs_stat.txt", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.gz_stat"), emit: quality_stats
  tuple val(sampleID), file("*_fastqs_stat.txt"), emit: fastqs_stats
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:
  
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
  cp *.gz_stat ${sampleID}_fastqs_stat.txt
  """
}
