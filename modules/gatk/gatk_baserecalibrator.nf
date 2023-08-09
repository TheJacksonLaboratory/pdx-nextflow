process GATK_BASERECALIBRATOR {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container 'quay.io/jaxpdx/gatk_3.4_java_1.7_samtools:latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.grp", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file("*.grp"), emit: grp

  script:
  
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
