process GATK_FILTERMUTECTCALLS {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '05:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/gatk-4.0.5.1_java_1.8_htslib_tabix.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf), file(tbi)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- Variant Calling FilterMutectCalls Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar \
  FilterMutectCalls \
  --variant ${sampleID}_intermed.vcf.gz \
  --output ${sampleID}_final_mutect_snp_indel_filtered.vcf \
  --unique-alt-read-count 5

  """
}
