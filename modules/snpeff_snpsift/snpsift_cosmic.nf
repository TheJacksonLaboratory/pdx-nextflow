process SNPSIFT_COSMIC {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container '/pdx/pdx_resource_service/elion/containers/java_perl_vcftools_python_2_snpEff_4_3.sif'


  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("${sampleID}_all_genes_variants_microindels_cosmicannotation.vcf"), emit: vcf

  script:

  """
  cat ${vcf} | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > ${sampleID}_onePerLine.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar annotate -id ${params.Cosmic_newer} ${sampleID}_onePerLine.vcf > ${sampleID}_all_genes_variants_microindels_cosmicannotation.vcf
  
  """
}