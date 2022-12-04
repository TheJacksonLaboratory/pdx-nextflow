process SUMMARY_STATS {
    tag "$sampleID"

    cpus = 1
    memory = 5.GB
    time = '01:00:00'

    container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/bedtools_2.27.1_python_2.7.3_java_1.8_GATK_3.4_samtools_1.3.1.sif'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'summary_stats' }", mode:'copy'

    input:
    tuple val(sampleID), file(filter_stat), file(duplicate_metrics), file(cov_metrics)

    output:
    file("${sampleID}_summary_stats.txt")

    script:
    """

    python  ${projectDir}/bin/wes/aggregate_stats_updated.py ${sampleID}_summary_stats.txt ${filter_stat} ${duplicate_metrics} ${cov_metrics}

    """
}