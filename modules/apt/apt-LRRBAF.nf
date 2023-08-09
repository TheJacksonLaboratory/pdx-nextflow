process LRRBAF {

    tag "$sampleID"

    cpus 1
    memory 6.GB
    time 8.h
    errorStrategy 'finish'

    container 'quay.io/jaxpdx/apt2.11.3_python2.7.11:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'apt' }", pattern: "*.{log}", mode:'copy'

    input:
        tuple val(sampleID), path(cel_list_file), path(cel_file)
        // Note: that the cel_file path is included here only to bring the cel file into the working directory of this step. 
        // The relative path is included in the `cel_list_file` and so it must be in the relative directory at runtime of this step.

    output:
        tuple val(sampleID), file("lrr_baf.txt"), emit: lrrbaf
        file("*log") //this might not be needed. Output was not originally being saved. 

    script:
    """
    /apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-probeset-summarize --cdf-file ${params.snp6chip} --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch ${params.hapmap_norm_target} --out-dir .  --cel-files ${cel_list_file}

    ${projectDir}/bin/cnv/normalize_affy_geno_cluster.pl ${params.genoclust} quant-norm.pm-only.med-polish.expr.summary.txt -locfile ${params.gw6_pfb_file} -out lrr_baf.txt

    """
}
//    rm -rf quant-norm.pm-only.med-polish.expr.summary.txt quant-norm.pm-only.med-polish.expr.report.txt 
//    This was originally included, but is not needed in the nextflow context. 
