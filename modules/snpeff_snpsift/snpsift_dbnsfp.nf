process SNPSIFT_DBNSFP{
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container '/pdx/pdx_resource_service/elion/containers/java_perl_vcftools_python_2_snpEff_4_3.sif'

  input:
  tuple val(sampleID), file(vcf)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  if (indel_snp == 'INDEL'){
    output_suffix = 'INDEL_snpsift_dbNSFPanno.vcf'
  }
  if (indel_snp =='SNP'){
    output_suffix = 'SNP_snpsift_dbNSFPanno.vcf'
  }
  if (indel_snp == 'BOTH'){
    output_suffix = 'all_genes_variants_microindels_snpEff_snpSift.vcf'
  }  
  if (indel_snp == 'UGPINDEL'){
    output_suffix = 'all_genes_variants_snpeff_snpsift.vcf'
  }
  """
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar dbnsfp -v \
  -db ${params.dbNSFP} -noDownload \
  -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AC,ExAC_AF \
  ${vcf} > ${sampleID}_${output_suffix}
  """
}