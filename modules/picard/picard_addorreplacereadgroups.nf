process PICARD_ADDORREPLACEREADGROUPS {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container '/pdx/pdx_resource_service/elion/containers/java_samtools_python_R_picard.sif'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(read_groups), file(bam)


  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar AddOrReplaceReadGroups \
  INPUT=${bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_groups.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_groups) \
  CREATE_INDEX=true
  """
}
