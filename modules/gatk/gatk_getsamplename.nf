process GATK_GETSAMPLENAME {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '05:00:00'

  container '/pdx/pdx_resource_service/elion/containers/gatk-4.0.5.1_java_1.8_htslib_tabix.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_tumor_SN.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file("*_tumor_SN.txt") 

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  java -Xmx${my_mem}G -jar /gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar  \
  GetSampleName \
  -I ${bam} \
  -O ${sampleID}_tumor_SN.txt
  """
}
