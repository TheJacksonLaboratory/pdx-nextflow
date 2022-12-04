process TMB_SCORE_CTP {
  tag "$sampleID"

  cpus 1
  memory 25.GB
  time '01:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/bedtools_2.27.1_java_1.8_snpeff_4.3_R.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'tmb' }", mode:'copy'

  input:
  tuple val(sampleID), file(variant_vcf), file(pindel_vcf)

  output:
  file "*.tab"
  file "*score"
  tuple val(sampleID), file("*HM.tab")
  tuple val(sampleID), file("*score")

  shell:
  
  '''

  cat !{variant_vcf} !{pindel_vcf} > !{sampleID}_variant_pindel.vcf

  bedtools intersect -header -v -a !{sampleID}_variant_pindel.vcf -b !{params.rec_var}/*bed | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "strandBias" --rmFilter "PASS" 'FS > 60' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowMQ" --rmFilter "PASS" 'MQ < 40' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowMQRankSum" --rmFilter "PASS" 'MQRankSum < -12.5' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowReadPosRankSum" --rmFilter "PASS" 'ReadPosRankSum < -8' > !{sampleID}_nofalsePositives.tmp2.vcf

  java -jar /snpEff_v4_3/snpEff/snpEff.jar eff -v -dataDir !{params.snpEff_data} -lof -canon -hgvs hg38 -noStats !{sampleID}_nofalsePositives.tmp2.vcf > !{sampleID}_all_genes_variants_snpEff.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar dbnsfp -v -db !{params.dbNSFP} -noDownload -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AC,ExAC_AF !{sampleID}_all_genes_variants_snpEff.vcf > !{sampleID}_all_genes_variants_snpeff_snpsift.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar annotate -id !{params.Cosmic_older} !{sampleID}_all_genes_variants_snpeff_snpsift.vcf > !{sampleID}_all_genes_variants_cosmicannotation.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "PutativeGermline" --rmFilter "PASS" '(ALT_AF[ANY] >= 90 & dbNSFP_1000Gp3_AF[ANY] >= 0.0095) | ((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_1000Gp3_AF[ANY] >= 0.0095)| (ALT_AF[ANY] >= 90 & dbNSFP_ESP6500_AA_AF[ANY] >= 0.0095)| ((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ESP6500_AA_AF[ANY] >= 0.0095)|(ALT_AF[ANY] >= 90 & dbNSFP_ESP6500_EA_AF[ANY] >= 0.0095)|((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ESP6500_EA_AF[ANY] >= 0.0095)|(ALT_AF[ANY] >= 90 & dbNSFP_ExAC_AF[ANY] >= 0.0095)|((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ExAC_AF[ANY] >= 0.0095)' !{sampleID}_all_genes_variants_cosmicannotation.vcf > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag.vcf

  cat !{sampleID}_all_genes_variants_cosmicannotation_germlineflag.vcf | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar extractFields !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline.vcf CHROM POS REF ALT ID FILTER DP_HQ ALT_AF "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_ExAC_AC" "dbNSFP_ExAC_AF" "CALLER" > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline_final.vcf


  header="chr\\tstart\\tend\\tCTPlength\\tHM"

  cat !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline_final.vcf | grep "HIGH\\|MODERATE" | awk -F '\\t' '{ if($5 !~ "rs" && ($6==""||$6=="PASS"||$6==".")) print $1,$2-1,$2 }'| sort| uniq| tr ' ' '\\t' > !{sampleID}_count2 ; echo -e ${header} > !{sampleID}_HM.tab


  bedtools coverage -a !{params.bins_ctpcoverage} -b !{sampleID}_count2 | cut -f 1-5 >> !{sampleID}_HM.tab


  Rscript !{projectDir}/bin/wes/TMB_final_CTP.R !{sampleID}_HM.tab !{sampleID}_TMB.score

  '''
}