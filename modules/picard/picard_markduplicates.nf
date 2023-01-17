process PICARD_MARKDUPLICATES {
  tag "$sampleID"

  cpus 1
  memory 16.GB
  time '12:00:00'

  container '/pdx/pdx_resource_service/elion/containers/picard_2.8.1.sif'

  // save if mouse and wes or save if keep intermediate
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.bam", mode:'copy', enabled: params.workflow=='wes' ? true : params.keep_intermediate
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "${sampleID}_dup_metrics.dat", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_dedup.bam"), emit: dedup_bam
  tuple val(sampleID), file("*_dedup.bai"), emit: dedup_bai
  tuple val(sampleID), file("${sampleID}_dup_metrics.dat"), emit: dedup_metrics

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /picard.jar MarkDuplicates \
  I=${bam} \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.dat \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT
  """
}
