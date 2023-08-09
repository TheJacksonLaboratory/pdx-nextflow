process SNPEFF_ANNOTATE {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container 'quay.io/jaxpdx/java_perl_vcftools_python_2_snpeff_4_3:latest'


  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("${sampleID}_all_genes_variants_microindels_snpEff.vcf"), emit: vcf

  script:

  """
  
  java -jar /snpEff_v4_3/snpEff/snpEff.jar eff -c ${params.snpEff_config} -v -lof -dataDir ${params.hgvs_data} -onlyTr ${params.ensembl_transcript} -hgvs GRCh38.84 -noStats ${vcf} > ${sampleID}_all_genes_variants_microindels_snpEff.vcf
  
  """
}