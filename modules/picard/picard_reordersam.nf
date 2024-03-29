process PICARD_REORDERSAM {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/jaxpdx/java_samtools_python_r_picard:latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar ReorderSam \
  INPUT=${bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  REFERENCE=${params.ref_fa} \
  CREATE_INDEX=true
  """
}
