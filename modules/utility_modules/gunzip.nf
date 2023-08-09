process GUNZIP {

  tag "$sampleID"

  cpus 1  
  memory { 5.GB * task.attempt }
  time { 2.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container "quay.io/jaxpdx/python_2.7.3:latest"

  stageInMode 'copy'

  input:
  tuple val(sampleID), file(reads)

  output:
  tuple val(sampleID), file("*.{fastq,fq}"), emit: gunzip_fastq
  shell:

  '''
  if [[ !{reads[0]} =~ ".gz" ]];
  then
    gunzip !{reads[0]}
  else
    mv !{reads[0]} input_!{reads[0]}
  fi
  if [[ !{reads[1]} =~ ".gz" ]];
  then
    gunzip !{reads[1]}
  else
   mv !{reads[1]} input_!{reads[1]}
  fi
  '''
}