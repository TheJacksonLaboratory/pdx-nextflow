#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {param_log} from "${projectDir}/bin/log/cnv.nf"
include {helpMessage} from "${projectDir}/bin/help/cnv.nf"
include {RUN_START} from "${projectDir}/bin/shared/run_start"

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

// main workflow
workflow CNV {
    RUN_START()
}