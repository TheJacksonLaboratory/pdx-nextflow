process JOIN_ADJACENT_SNPS_AS {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container '/pdx/pdx_resource_service/elion/containers/perl_java_1.8_python_2.7.3_tabix_samtools_bcftools_htslib_snpEff_4_3.sif'


  input:
  tuple val(sampleID), file(bam), file(bai), file(filt_var_DP)

  output:
  tuple val(sampleID), file("*fixAdjSNP.vcf.gz"), emit: vcf
  tuple val(sampleID), file("*fixAdjSNP.vcf.gz.tbi"), emit: tbi

  script:
  
  """
  python ${projectDir}/bin/wes/joinAdjacentSNPs_AS.py -v ${filt_var_DP} \
  -o ${sampleID}_variant_fixAdjSNP.vcf \
  1 ${bam} ${params.fa2bit} 2>${sampleID}_AdjSNP.log.txt

  awk 'NF' ${sampleID}_variant_fixAdjSNP.vcf > temp

  mv temp ${sampleID}_variant_fixAdjSNP.vcf

  bgzip ${sampleID}_variant_fixAdjSNP.vcf

  tabix ${sampleID}_variant_fixAdjSNP.vcf.gz

  """
}
