process GATK_INDELREALIGNER{
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '08:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/GATK_3.4_java_1.7_samtools.sif'

  // save if mouse and wes or save if keep intermediate
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'gatk' }", pattern: "*.bam", mode:'copy', enabled: params.workflow=='wes' ? true : params.keep_intermediate

  input:
  tuple val(sampleID), file(bam), file(bai), file(intervals)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  log.info "----- GATK IndelRealigner Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /GenomeAnalysisTK.jar \
  -I ${bam} \
  -R ${params.ref_fa} \
  -T IndelRealigner \
  -targetIntervals ${intervals} \
  -o ${sampleID}_realigned.bam \
  -known ${params.gold_std_indels} \
  -known ${params.known_indels} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -L ${params.targets_gatk}
  """
}
