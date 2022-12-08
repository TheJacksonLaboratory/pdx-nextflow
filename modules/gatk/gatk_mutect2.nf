process GATK_MUTECT2 {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '05:00:00'

  container '/pdx/pdx_resource_service/elion/containers/gatk-4.0.5.1_htslib_tabix.sif'

  input:
  tuple val(sampleID), file(bam), file(bai), file(tumor)

  output:
  tuple val(sampleID), file("*.vcf.gz"), emit: vcf
  tuple val(sampleID), file("*.vcf.gz.tbi"), emit: tbi

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  tumorName=\$(cat ${tumor})

  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /gatk/gatk.jar \
  Mutect2 \
  -R ${params.ref_fa} \
  -I ${bam}  \
  -tumor \${tumorName} \
  --germline-resource ${params.exac_ref} \
  --af-of-alleles-not-in-resource 0.0000082364 \
  -O ${sampleID}_intermed.vcf \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --dont-use-soft-clipped-bases false \
  --genotype-germline-sites false \
  --sample-ploidy ${params.samp_ploidy} \
  -L ${params.targets_gatk} \
  --annotation QualByDepth \
  --annotation RMSMappingQuality \
  --annotation FisherStrand \
  --annotation MappingQualityRankSumTest \
  --annotation ReadPosRankSumTest \
  --min-base-quality-score 20 \
  --standard-min-confidence-threshold-for-calling 30

  bgzip ${sampleID}_intermed.vcf

  tabix ${sampleID}_intermed.vcf.gz

  """
}
