process BWA_MEM {
  tag "$sampleID"

  cpus 12
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container '/pdx/pdx_resource_service/elion/containers/bwakit_v0.7.15_cv1.sif'

  input:
  tuple val(sampleID), file(fq_reads), file(read_groups)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam

  script:
  
  if (params.read_type == "SE"){
    inputfq="${fq_reads[0]}"
    }
  if (params.read_type == "PE"){
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

  """
  rg=\$(cat $read_groups)
  run-bwamem -R \${rg} \
  -t $task.cpus -o ${sampleID} -H ${params.ref_fa_indices} $inputfq | sh
  """
}
