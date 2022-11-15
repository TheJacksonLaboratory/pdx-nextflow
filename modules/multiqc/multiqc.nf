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