process GATK_PRINTREADS {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container 'quay.io/jaxpdx/gatk_3.4_java_1.7_samtools:latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'gatk' }", pattern: "*.bam", mode:'copy', enabled: params.workflow=='wes' ? true : params.keep_intermediate

  input:
  tuple val(sampleID), file(bam), file(bai), file(grp)

  output:
  tuple val(sampleID), file("*realigned_BQSR.bam"), emit: bam
  tuple val(sampleID), file("*realigned_BQSR.bam.bai"), emit: bai

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  java "-Xmx${my_mem}G" -jar /GenomeAnalysisTK.jar \
  -T PrintReads \
  -I ${bam} \
  -R ${params.ref_fa} \
  -BQSR ${grp} \
  -o ${sampleID}_realigned_BQSR.bam \
  --disable_auto_index_creation_and_locking_when_reading_rods

  samtools index ${sampleID}_realigned_BQSR.bam
  """
}
