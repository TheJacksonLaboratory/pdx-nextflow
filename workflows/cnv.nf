#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {param_log} from "${projectDir}/bin/log/cnv.nf"
include {helpMessage} from "${projectDir}/bin/help/cnv.nf"
include {RUN_START} from "${projectDir}/bin/shared/run_start"
include {GET_MODEL_GENDER} from "${projectDir}/modules/apt/apt-get_model_gender"
include {LRRBAF} from "${projectDir}/modules/apt/apt-LRRBAF"
include {ASCAT} from "${projectDir}/modules/ascat/ascat"
//include {ASCAT_ANNOTATION} from "${projectDir}/modules/ascat/ascat-annotation"

// log params
param_log()

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
////// Required parameters \\\\\\
if ( ! params.celInput ) {
    exit 1, "Parameter ERROR: celInput ($params.celInput) must be the absolute path to CEL file."
}
if ( ! params.hapmap_dat ) {
  exit 1, "Parameter ERROR: Directory containing hapmap data from NCBI must be specified."
}
if ( ! params.hapmap_fm ) {
  exit 1, "Parameter ERROR: file containing hapmap data for females must be specified."
}
if ( ! params.hapmap_m ) {
  exit 1, "Parameter ERROR: file containing hapmap data for males must be specified."
}
if ( ! params.genoclust ) {
  exit 1, "Parameter ERROR: genotype clustering file must be specified."
}
if ( ! params.snp6chip ) {
  exit 1, "Parameter ERROR: file containing description of the probe sets on the chip must be specified."
}
if ( ! params.snp6chip_birdseed_mod ) {
  exit 1, "Parameter ERROR: file containing birdseed models must be specified."
}
if ( ! params.snp6chip_specsnps ) {
  exit 1, "Parameter ERROR: file containing chromosme X (non-pseudo-autosomal), chromosome Y, and mitochondrial SNPs must be specified."
}
if ( ! params.hapmap_norm_target ) {
  exit 1, "Parameter ERROR: file containing hapmap sample data used for normalization must be specified."
}
if ( ! params.gw6_pfb_file ) {
  exit 1, "Parameter ERROR: file containing annotated marker positions must be specified."
}
if ( ! params.SNPpos ) {
  exit 1, "Parameter ERROR: file containing annotated marker positions must be specified."
}
if ( ! params.GC ) {
  exit 1, "Parameter ERROR: file containing GC% for annotated marker positions must be specified."
}
if ( ! params.chr_arm ) {
  exit 1, "Parameter ERROR: file containing positions of chromosome arms must be specified."}
if ( ! params.exp_mart_genes ) {
  exit 1, "Parameter ERROR: file containing genes exported to PDX data mart database must be specified."
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def celRE = ~/.CEL$/

celin = params.celInput
celin2 = file(celin)

// Need to specify absolute path to CEL file because relative paths cause symlink breakage
cel = celin2.toAbsolutePath().toString()

def sampleID = ( celin2.name - celRE )


//~~~~~~~~~~ Initial Channel of SampleID and CEL Data ~~~~
Channel.of( sampleID, cel )
       .toList()
       .set { sample_CEL_ch }

// main workflow
workflow CNV {

    // Remove `pipeline_complete.txt` from prior run, if this is a 'resume' or sample re-run. 
    // This file is used in 'on.complete' and in JAX PDX loader
    run_check = file("${params.pubdir}/pipeline_complete.txt")
    run_check.delete()
    
    RUN_START()
    
    GET_MODEL_GENDER(sample_CEL_ch)

    lrrbaf_input = GET_MODEL_GENDER.out.cel_list.join.sample_CEL_ch

    LRRBAF(lrrbaf_input)

    ascat_input = LRRBAF.out.lrrbaf.join.GET_MODEL_GENDER.out.gender

    ASCAT(ascat_input)

    ascat_annotation_input = ASCAT.out.raw_seg.join.ASCAT.out.ploidy.join.GET_MODEL_GENDER.out.gender

    ASCAT_ANNOTATION(ascat_annotation_input)

    // rm *.RData *.BAF.txt *.LogR.txt lrr_baf1.txt *.PCFed.txt

}

workflow.onComplete {
    if (workflow.success && params.preserve_work == "no") {
    workflow.workDir.deleteDir()
    log.info "Cleaned Work Directory"
    } else {
    log.info "Keeping Work Directory"
    }
    if (workflow.success) {
    log.info "Pipeline completed successfully"
    run_check = file("${params.pubdir}/pipeline_running.txt")
    run_check.renameTo("${params.pubdir}/pipeline_complete.txt")
    } else {
    log.info "Pipeline completed with errors"
    }
}
