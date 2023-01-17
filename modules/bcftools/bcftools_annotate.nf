process BCF_ANNOTATE {
  tag "$sampleID"

  cpus 1
  memory 6.GB
  time '06:00:00'

  container '/pdx/pdx_resource_service/elion/containers/perl_java_1.8_python_2.7.3_tabix_samtools_bcftools_htslib_snpEff_4_3.sif'


  input:
  tuple val(sampleID), file(vcf), file(tbi)

  output:
  tuple val(sampleID), file("*.noIds.vcf.gz"), emit: vcf

  script:
  
  """
  bcftools annotate \
  --output ${sampleID}_variant_fixAdjSNP.noIds.vcf.gz \
  --output-type z \
  --remove ID ${vcf}
  """

}
