process ASCAT_ANNOTATION {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time 4.h
    errorStrategy 'finish'

    container '/projects/omics_share/.pdx/pdx_resource_service/elion/containers/devtools_ASCAT_R.sif'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'ascat' }", pattern: "*segments*.txt", mode:'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'ascat' }", pattern: "*png", mode:'copy'
    
    input:
        tuple val(sampleID), path(raw_segs), path(sample_ploidy), path(sample_gender)
 
    output:
        tuple val(sampleID), file("*segments_raw.extend.ensgene_cnvbreak.txt"), emit: ascast_annot
        file("*segments*txt")
        file("*png")

    shell:
    '''

    !{projectDir}/bin/cnv/segment_raw_extend.pl !{raw_segs} !{sample_ploidy} !{params.chr_arm} !{sample_gender}

    !{projectDir}/bin/cnv/ensemblegenes_cnv_break.pl !{sampleID}.segments_raw.extend.txt !{params.exp_mart_genes}

    !{projectDir}/bin/cnv/get_msp.py !{sampleID}.segments_raw.extend.txt > msp.txt

    msp=$(cat msp.txt)

    R CMD BATCH --slave "--args !{sampleID}.segments_raw.extend.txt $msp ./ " !{projectDir}/bin/cnv/seg_plot.R

    '''

}

//     rm -rf tmp.txt  
//     This was originally included, but is not needed in the nextflow context. 