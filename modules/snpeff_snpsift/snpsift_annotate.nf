process SNPSIFT_ANNOTATE{
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/perl_java_1.8_python_2.7.3_tabix_snpeff_4_3.sif'


  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf*"), emit: vcf

  script:
  log.info "----- snpSift Annotate Running on: ${sampleID} -----"

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  
  if (vcf =~ 'mutect'){
    options = 'mutect_snp_indel_filtered.vcf.additionalfilters.tmp.vcf'
  }
  if (vcf =~ 'SNP'){
    options = 'variant_fixAdjSNP.vcf'
  }
  if (vcf =~ 'microIndels'){
    options = 'microIndels.DPfiltered1.vcf.additionalfilters.tmp.vcf'
  }

  """
  java -Xmx${my_mem}G -jar /snpEff_v4_3/snpEff/SnpSift.jar \
  annotate \
  -id ${params.dbsnp} ${vcf} > ${sampleID}_${options}
  """
}
