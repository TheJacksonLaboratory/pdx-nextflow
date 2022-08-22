process CLASSIFICATION_A {

  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '24:00:00'

  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/xenome_1.0.1.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'xenome' }", pattern: "*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(trimmed)

  output:
  tuple val(sampleID), file("human*{1,2}.fastq"), emit: xenome_fastq
  tuple val(sampleID), file("*.txt"), emit: xenome_stats

  script:
  log.info "----- Xenome running on: ${sampleID} -----"
  """
  /xenome-1.0.1-r/xenome classify -T 12 -P ${params.ref_prefix} --pairs --host-name mouse --graft-name human -i ${trimmed[0]} -i ${trimmed[1]} > ${sampleID}_xenome_stats.txt

  rm -rf *both*fastq* *mouse*fastq* *neither*fastq* *ambiguous*fastq*
  """
}
