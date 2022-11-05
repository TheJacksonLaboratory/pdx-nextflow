process SNPSIFT_MICROINDELS_CTP {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/java_perl_vcftools_python_2_snpEff_4_3.sif'


  input:
  tuple val(sampleID), file(variants_vcf), file(microindels_vcf)

  output:
  tuple val(sampleID), file("${sampleID}_all_genes_variants_microindels_filtered.vcf"), emit: vcf

  script:

  """
  
  cat ${variants_vcf} ${microindels_vcf} > ${sampleID}.merged.unsorted.vcf
  vcf-sort ${sampleID}.merged.unsorted.vcf > ${sampleID}_all_genes_variants_microindels.vcf
  
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' ${sampleID}_all_genes_variants_microindels.vcf | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "PASS" --rmFilter "clustered_events" "( FILTER == 'clustered_events' )" | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;minDP' )" | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk' )" | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk;minDP' )" > ${sampleID}_all_genes_variants_microindels_filtered.vcf

  """
}