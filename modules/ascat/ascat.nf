process ASCAT {

    tag "$sampleID"

    cpus 1
    memory 2.GB
    time 4.h
    errorStrategy 'finish'

    container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/devtools_ASCAT_R.sif'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'ascat' }", pattern: "*.{txt,png,Rout}", mode:'copy'

    input:
        tuple val(sampleID), path(lrr_baf), path(gender)

    output:
        tuple val(sampleID), file("*segments_raw.txt"), emit: raw_seg
        tuple val(sampleID), file("*ploidy.txt"), emit: ploidy
        file("*txt")
        file("*png")
        file("*Rout")
        // file("*RData") this could be saved as an intermediate. It was not ultimately saved in the pipeline. 

    shell:
    '''

    cut -f1-5 !{lrr_baf} > lrr_baf1.txt

    R CMD BATCH --slave "--args !{params.SNPpos} !{params.GC}"  !{projectDir}/bin/cnv/lrrbaf_ascat_tumor.R

    '''
}

