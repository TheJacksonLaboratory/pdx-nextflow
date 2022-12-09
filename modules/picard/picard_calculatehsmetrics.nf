process PICARD_CALCULATEHSMETRICS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/picard-1.95_python_2_7_3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file("${sampleID}_CoverageMetrics.txt"), emit: hsmetrics

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Djava.io.tmpdir=$TMPDIR -jar -Xmx${my_mem}G /picard-tools-1.95/CalculateHsMetrics.jar \
  TARGET_INTERVALS=${params.targets_picard} \
  BAIT_INTERVALS=${params.targets_exons_picard} \
  REFERENCE_SEQUENCE=${params.ref_fa} \
  INPUT=${bam} \
  OUTPUT=${sampleID}_CoverageMetrics.txt \
  VALIDATION_STRINGENCY=LENIENT

  """
}
