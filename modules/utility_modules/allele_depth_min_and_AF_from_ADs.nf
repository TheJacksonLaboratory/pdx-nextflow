process ALLELE_DEPTH_MIN_AND_AF_FROM_ADS {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container "/pdx/pdx_resource_service/elion/containers/python_2.7.3.sif"

  input:
  tuple val(sampleID), file(mutect2_filtered)

  output:
  tuple val(sampleID), file("*DPfiltered.tmp.vcf"), emit: vcf

  script:
  
  """
  python ${projectDir}/bin/exome/allele_depth_min_and_AF_from_ADs.py ${mutect2_filtered} ${sampleID}_mutect_snp_indel_filtered.vcf.DPfiltered.tmp.vcf ${params.minDP}
  """
}
