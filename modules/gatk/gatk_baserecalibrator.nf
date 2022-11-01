process GATK_BASERECALIBRATOR {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/GATK_3.4_java_1.7_samtools.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'gatk' }", pattern: "*.grp", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*.grp"), emit: grp

  script:
  log.info "----- GATK BaseRecalibrator Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  java -Djava.io.tmpdir=$TMPDIR "-Xmx${my_mem}G" -jar /GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -I ${bam} \
  -R ${params.ref_fa} \
  -knownSites ${params.dbsnp} \
  -knownSites ${params.gold_std_indels} \
  -knownSites ${params.known_indels} \
  -o ${sampleID}_recal_data.grp \
  --disable_auto_index_creation_and_locking_when_reading_rods
  """
}
