process COMPRESS_INDEX_VCF {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container '/pdx/pdx_resource_service/elion/containers/quay.io-biocontainers-samtools-1.14--hb421002_0.img'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", mode:'copy'

  input:
    tuple val(sampleID), file(vcf)

  output:
      tuple val(sampleID), file("*.vcf.gz"), emit: vcf
      tuple val(sampleID), file("*.tbi"), emit: tbi

  script:
    """
    bgzip ${sampleID}_intermed.vcf

    tabix ${sampleID}_intermed.vcf.gz
    """
}
