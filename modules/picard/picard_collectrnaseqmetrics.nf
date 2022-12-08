process PICARD_COLLECTRNASEQMETRICS {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '03:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/java_samtools_python_R_picard.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.bam", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.pdf", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*picard_aln_metrics.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics
  tuple val(sampleID), file("*.pdf"), optional: true

  script:
  
    if (params.read_prep == "reverse_stranded") {
    strand_setting = "SECOND_READ_TRANSCRIPTION_STRAND"
  }

  if (params.read_prep == "forward_stranded") {
    strand_setting = "FIRST_READ_TRANSCRIPTION_STRAND"
  }

  if (params.read_prep == "non_stranded") {
    strand_setting = "NONE"
  }

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar CollectRnaSeqMetrics \
  I=${bam} \
  O=${sampleID}_picard_aln_metrics.txt \
  REF_FLAT=${params.ref_flat} \
  RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
  STRAND=${strand_setting} \
  CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
  """
}
