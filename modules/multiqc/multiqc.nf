process MULTIQC {

    container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/quay.io-biocontainers-multiqc-1.12--pyhdfd78af_0.img'
    
    publishDir "${params.pubdir}/multiqc", pattern: "*multiqc_report.html", mode:'copy'
    publishDir "${params.pubdir}/multiqc", pattern: "*_data", mode:'copy'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data" , emit: data
    path "*_plots" , optional:true, emit: plots

    script:

    """
    multiqc .
    """

}


  process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:

    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])

    file (fusions_mq) from summary_fusions_mq.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    when: !params.debug

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}