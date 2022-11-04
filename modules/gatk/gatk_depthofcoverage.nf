process GATK_DEPTHOFCOVERAGE {

  tag "$sampleID"

  cpus 2
  memory 15.GB
  time '05:00:00'


  container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/bedtools_2.27.1_python_2.7.3_java_1.8_GATK_3.4_samtools_1.3.1.sif'


  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)
  val(L)

  output:
  tuple val(sampleID), file("*_gatk_temp*.txt"), emit: txt

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  id = L =~ /probes/ ?  "1" : "4"
  """
  java -Djava.io.tmpdir=$TMPDIR "-Xmx${my_mem}g" -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp${id}.txt \
  -I ${bam} \
  -L  ${L} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS
  """
}
