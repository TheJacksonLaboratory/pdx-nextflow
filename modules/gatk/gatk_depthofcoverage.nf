process GATK_DEPTHOFCOVERAGE {

  tag "$sampleID"

  cpus 2
  memory 15.GB
  time '05:00:00'


  container '/pdx/pdx_resource_service/elion/containers/bedtools_2.27.1_python_2.7.3_java_1.8_GATK_3.4_samtools_1.3.1.sif'

  input:
  tuple val(sampleID), file(bam), file(bai)
  val(L)

  output:
  tuple val(sampleID), file("${sampleID}_gatk_coverage_formatted.txt"), emit: txt
  
  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  id = L =~ /PROBES/ ?  "1" : "4"

  if (params.workflow == 'rnaseq')
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
  
  cp ${sampleID}_gatk_temp${id}.txt ${sampleID}_gatk_coverage_formatted.txt
  """

  else if (params.workflow == 'ctp')
  """
  java -Djava.io.tmpdir=$TMPDIR "-Xmx${my_mem}g" -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp1.txt \
  -I ${bam} \
  -L ${L} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable

  bash ${projectDir}/bin/exome/gatk_formatter.sh ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt ${L}
  cp ${sampleID}_gatk_temp3.txt ${sampleID}_gatk_coverage_formatted.txt
  """

  else if (params.workflow == 'wes')
  """
  java -Djava.io.tmpdir=$TMPDIR "-Xmx${my_mem}g" -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp1.txt \
  -I ${bam} \
  -L ${L} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable

  bash ${projectDir}/bin/exome/gatk_formatter.sh ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt ${L}
  cp ${sampleID}_gatk_temp3.txt ${sampleID}_gatk_coverage_formatted.txt
  """
  else
  error "Invalid workflow: $params.workflow"

}
