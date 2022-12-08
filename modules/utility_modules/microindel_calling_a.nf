process MICROINDEL_CALLING_A {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container '/pdx/pdx_resource_service/elion/containers/pindel_v0.2.5.sif'


  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  
  """
  echo -e "${bam}\t350\t${sampleID}" > ${sampleID}_pindel_config.txt

  pindel \
  --fasta ${params.ref_fa} \
  --config-file ${sampleID}_pindel_config.txt \
  -o ${sampleID}

  cat  ${sampleID}_D  ${sampleID}_SI >  ${sampleID}_DSI

  /pindel-0.2.5/pindel2vcf \
  -p  ${sampleID}_DSI \
  -r ${params.ref_fa} \
  -R hg38 \
  -d 20150925 \
  --max_size 50 \
  --vcf ${sampleID}_microIndels_all.vcf \
  -G \
  --het_cutoff 0.05

  """
}
