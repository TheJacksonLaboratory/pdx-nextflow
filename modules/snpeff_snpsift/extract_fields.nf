process EXTRACT_FIELDS {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container '/pdx/pdx_resource_service/elion/containers/java_perl_vcftools_python_2_snpEff_4_3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'annotations' }", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)

  output:
  file("${sampleID}_variants_microIndels.DPfiltered.Annotated.tab")
  file("${sampleID}_variants_microIndels.Hardfiltered.Annotated.txt")


  shell:

  '''
  
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar extractFields !{vcf} CHROM POS REF ALT ID FILTER DP_HQ ALT_AF "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_ExAC_AC" "dbNSFP_ExAC_AF" "CALLER" > !{sampleID}_variants_microIndels.DPfiltered.Annotated.tab

  python !{projectDir}/bin/exome/clean_intergenic_region_gene_names.py -f !{sampleID}_variants_microIndels.DPfiltered.Annotated.tab -d

  cat !{sampleID}_variants_microIndels.DPfiltered.Annotated.tab | awk -F '\t' 'BEGIN {OFS="\t"} $6 == "FILTER" || $6 == "PASS" || $6 == "" || $6 == "."' > !{sampleID}_variants_microIndels.Hardfiltered.Annotated.txt
  
  '''
}