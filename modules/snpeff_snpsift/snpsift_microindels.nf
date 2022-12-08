process SNPSIFT_MICROINDELS {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container '/pdx/pdx_resource_service/elion/containers/java_perl_vcftools_python_2_snpEff_4_3.sif'

  input:
    tuple val(sampleID), file(variants_vcf), file(microindels_vcf)

  output:
    tuple val(sampleID), file("${sampleID}_all_genes_variants_microindels_filtered.vcf"), emit: vcf

  script:
  if (params.workflow == "ctp")
    """

    cat ${variants_vcf} ${microindels_vcf} > ${sampleID}.merged.unsorted.vcf
    vcf-sort ${sampleID}.merged.unsorted.vcf > ${sampleID}_all_genes_variants_microindels.vcf
    
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' ${sampleID}_all_genes_variants_microindels.vcf | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "PASS" --rmFilter "clustered_events" "( FILTER == 'clustered_events' )" | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;minDP' )" | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk' )" | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk;minDP' )" > ${sampleID}_all_genes_variants_microindels_filtered.vcf

    """
  else if (params.workflow == "wes")
    """

    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowQD" 'QD < 2.0'  ${variants_vcf} | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter  "strandbias" 'FS > 60.0' | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter  "lowMQ" 'MQ < 40'  | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter  "lowMQRankSum" 'MQRankSum < -12.5'  | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter  "lowReadPosRankSum" 'ReadPosRankSum < -8' > ${sampleID}_variants.vcf
    
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowQD" 'QD < 2.0' < ${microindels_vcf} | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter  "strandbias" 'FS > 200.0' | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter  "lowReadPosRankSum" 'ReadPosRankSum < -20' > ${sampleID}_microindels.vcf 


    cat ${sampleID}_variants.vcf ${sampleID}_microindels.vcf  > ${sampleID}.merged.unsorted.vcf

    vcf-sort ${sampleID}.merged.unsorted.vcf > ${sampleID}_all_genes_variants_microindels.vcf
    
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' ${sampleID}_all_genes_variants_microindels.vcf | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "PASS" --rmFilter "clustered_events" "( FILTER == 'clustered_events' )" | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;minDP' )" | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk' )" | \
    java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk;minDP' )" > ${sampleID}_all_genes_variants_microindels_filtered.vcf

    """
}