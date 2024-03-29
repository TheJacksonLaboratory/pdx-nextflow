process RNA_SUMMARY_STATS {
    tag "$sampleID"

    cpus = 1
    time = '00:15:00'

    container 'quay.io/jaxpdx/perl:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'summary_stats' }", pattern: "*stats.txt", mode:'copy'

    input:
    tuple val(sampleID), file(rsem_stats), file(quality_stats), file(xenome_stats), file(picard_metrics)

    output:
    tuple val(sampleID), file("*.txt")

    script:
    
    if (params.read_type == "PE")

      """
      perl ${projectDir}/bin/rnaseq/summary_QC_metrics.pl \
      ${quality_stats} \
      ${xenome_stats} \
      ${rsem_stats} \
      ${picard_metrics} > ${sampleID}_summary_stats.txt
      """

    else if (params.read_type == "SE")

      """
      perl ${projectDir}/bin/rnaseq/summary_QC_metrics.pl \
      ${quality_stats} \
      ${xenome_stats} \
      ${rsem_stats} \
      ${picard_metrics}  > ${sampleID}_summary_stats.txt
      """
}
