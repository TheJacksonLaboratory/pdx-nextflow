process BCF_ANNOTATE {
  tag "$sampleID"

  cpus 1
  memory 6.GB
  time '06:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/perl_java_1.8_python_2.7.3_tabix_samtools_bcftools_htslib_snpEff_4_3.sif'


  input:
  tuple val(sampleID), file(vcf), file(tbi)

  output:
  tuple val(sampleID), file("*.noIds.vcf.gz"), emit: vcf

  script:
  log.info "----- BCFTools Annotate Running on: ${sampleID} -----"

  """
  bcftools annotate \
  --output ${sampleID}_variant_fixAdjSNP.noIds.vcf.gz \
  --output-type z \
  --remove ID ${vcf}
  """

}
