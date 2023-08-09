process SUMMARY_STATS {
    tag "$sampleID"

    cpus = 1
    memory = 5.GB
    time = '01:00:00'

    container 'quay.io/jaxpdx/bedtools_2.27.1_python_2_7_3_java_1.8_gatk_3.4_samtools_1.3.1:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'summary_stats' }", mode:'copy'

    input:
    tuple val(sampleID), file(filter_stat), file(duplicate_metrics), file(cov_metrics)

    output:
    file("${sampleID}_summary_stats.txt")

    script:
    """

    python  ${projectDir}/bin/exome/aggregate_stats_updated.py ${sampleID}_summary_stats.txt ${filter_stat} ${duplicate_metrics} ${cov_metrics}

    """
}