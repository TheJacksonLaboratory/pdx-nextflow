#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show()



//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Configurable variable parameters specific to individual runs:
    params.celInput               = null // absolute path to input CEL.
    params.outdir                 = null // directory where output will be stored.
    params.tmpdir                 = "$TMPDIR"
    params.hapmap_dat             = null // directory containing hapmap data from NCBI
    params.hapmap_fm              = null // file containing hapmap data for females
    params.hapmap_m               = null // file containing hapmap data for males
    params.genoclust              = null // genotype clustering file
    params.snp6chip               = null // file containing description of the probe sets on the chip
    params.snp6chip_birdseed_mod  = null // file containing birdseed models
    params.snp6chip_specsnps      = null // file containing chromosme X (non-pseudo-autosomal), chromosome Y, and mitochondrial SNPs
    params.hapmap_norm_target     = null // file containing hapmap sample data used for normalization
    params.gw6_pfb_file           = null // file containing annotated marker positions
    params.SNPpos                 = null // file containing annotated marker positions
    params.GC                     = null // file containing GC% for annotated marker positions
    params.chr_arm                = null // file containing positions of chromosome arms
    params.exp_mart_genes         = null // file containing genes exported to PDX data mart database
}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow -c path/to/params.cfg run path/to/pipeline.nf -profile slurm, singularity
            (The params.cfg file needs to have the following mandatory parameters
             OR they need to specified on the command line.)

    Mandatory:
        --celInput                absolute path to input CEL.
        --outdir                  directory where output will be stored.
        --tmpdir                  = "$TMPDIR"
        --hapmap_dat              directory containing hapmap data from NCBI
        --hapmap_fm               file containing hapmap data for females
        --hapmap_m                file containing hapmap data for males
        --genoclust               genotype clustering file
        --snp6chip                file containing description of the probe sets on the chip
        --snp6chip_birdseed_mod   file containing birdseed models
        --snp6chip_specsnps       file containing chromosme X (non-pseudo-autosomal), chromosome Y, and mitochondrial SNPs
        --hapmap_norm_target      file containing hapmap sample data used for normalization
        --gw6_pfb_file            file containing annotated marker positions
        --SNPpos                  file containing annotated marker positions
        --GC                      file containing GC% for annotated marker positions
        --chr_arm                 file containing positions of chromosome arms
        --exp_mart_genes          file containing genes exported to PDX data mart database


    Optional:
        --tmpdir                Path to hold temporary files for software execution.
        -name                   Name for the pipeline run. If not specified Nextflow will
                                automatically generate a random mnemonic.

    """.stripIndent()
}

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
if ( ! params.outdir ) {
    exit 1, "Parameter ERROR: Output directory parameter must be specified."
}
if ( ! file(params.outdir).exists() ) {
    exit 1, "Parameter ERROR: Missing output directory ($params.outdir) is not found: check if path is correct."
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
  exit 1, "Parameter ERROR: file containing positions of chromosome arms must be specified."
}
if ( ! params.exp_mart_genes ) {
  exit 1, "Parameter ERROR: file containing genes exported to PDX data mart database must be specified."
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
def summary = [:]
summary['Pipeline']         = workflow.manifest.name
summary['Description']      = workflow.manifest.description
if(workflow.revision) {
    summary['Pipeline Release'] = workflow.revision
}
summary['Run Name']         = workflow.runName
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Config Files']     = workflow.configFiles
summary['Command Line']     = workflow.commandLine
summary['Nextflow Info']    = "v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Workflow dir']     = workflow.projectDir

// Pipeline Params:
summary['Parameters......']                                   = ''
summary['.  Output dir']                                      = params.outdir
summary['.  celInput']                                        = params.celInput
summary['.  hapmap data']                                     = params.hapmap_dat
summary['.  hapmap female']                                   = params.hapmap_fm
summary['.  hapmap male']                                     = params.hapmap_m
summary['.  genotype clustering']                             = params.genoclust
summary['.  probe sets on chip']                              = params.snp6chip
summary['.  birdseed models']                                 = params.snp6chip_birdseed_mod
summary['.  chrX (non-pseudo-autosomal), chrY, and MT SNPs']  = params.snp6chip_specsnps
summary['.  hapmap sample data for normalization']            = params.hapmap_norm_target
summary['.  annotated marker positions']                      = params.gw6_pfb_file
summary['.  annotated marker positions']                      = params.SNPpos
summary['.  GC% for annotated marker positions']              = params.GC
summary['.  chromosome arm positions']                        = params.chr_arm
summary['.  genes exported to PDX data mart']                 = params.exp_mart_genes
summary['.  TMP dir']                                         = params.tmpdir
summary['Run Start Time']                                     = workflow.start

//print summary header:
println "$summary"


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def tmpdir    = params.tmpdir


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


//~~~~~~~~~~~~~~~~ Primary publishDir == sample_outdir ~~~~~
def sample_outdir = "${params.outdir}/${sampleID}"
def sample_outputdir = "${params.outputdir}/${sampleID}"



// Step 1: Get Model Gender
process hCNV_get_model_gender {
  tag "sampleID"
  label 'med_mem'
  label 'apt2_11_3_python2_7_11'

  publishDir "${sample_outputdir}_tmp", pattern: "*.txt", mode: 'copy'
  publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'

  input:
  tuple sampleID, cel_file from sample_CEL_ch

  output:
  tuple sampleID, file("listfile1") into cel_list, dummy_cel_list
  tuple sampleID, file("gender.txt") into gender1, gender2, dummy_gender
  file("*txt")

  script:
  """

  sname=\$(echo ${sampleID} | cut -d"_" -f1)

  python ${params.model_gender} \${sname} > gender_elims.txt

  echo "cel_files" >> listfile

  echo $cel >> listfile

  awk -v a=${params.hapmap_dat} '{ if (NR>1) print a"/"\$2}' ${params.hapmap_fm} >> listfile

  awk -v a=${params.hapmap_dat} '{ if (NR>1) print a"/"\$2}' ${params.hapmap_m} >> listfile

  /apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-probeset-genotype -c ${params.snp6chip} -a birdseed --read-models-birdseed ${params.snp6chip_birdseed_mod} --special-snps ${params.snp6chip_specsnps} --out-dir ${sample_outputdir}_tmp  --cel-files listfile

  cat ${sample_outputdir}_tmp/birdseed.report.txt | grep -v "#" | awk 'NR==2' | cut -f2 > gender_birdseed.txt

  gender=\$(echo gender_elims.txt)

  if [ "X\$gender" = "Xunknown" -o "X\$gender" = "Xunspecified" ]; then gender=\$(echo gender_birdseed.txt); cp gender_birdseed.txt gender.txt; else cp gender_elims.txt gender.txt; fi

  echo "cel_files" > listfile1

  echo $cel >> listfile1

  if [ "X\$gender" = "Xfemale" -o "X\$gender" = "Xunknown" ]; then awk -v a=${params.hapmap_dat} '{if (NR>1) print a"/"\$2}' ${params.hapmap_fm} >> listfile1; elif [ "X\$gender" = "Xmale" ]; then awk -v a={params.hapmap_dat} '{if (NR>1) print a"/"\$2}' ${params.hapmap_m} >> listfile1; fi

  cat ${sample_outputdir}_tmp/birdseed.confidences.txt | grep -v "#" | cut -f1-2 > birdseed.confidences1.txt

  cat ${sample_outputdir}_tmp/birdseed.calls.txt | grep -v "#" | cut -f1-2 > birdseed.calls1.txt

  cat ${sample_outputdir}_tmp/birdseed.report.txt | grep -v "#" | head -2 > birdseed.report1.txt

  rm -rf listfile ${sample_outputdir}_tmp/birdseed.report.txt ${sample_outputdir}_tmp/birdseed.confidences.txt ${sample_outputdir}_tmp/birdseed.calls.txt

  """

  }

// Step 2: LRRBAF
process hCNV_LRRBAF {
  tag "sampleID"
  label 'med_mem'
  label 'apt_perl'


  publishDir "${sample_outputdir}_tmp", pattern: "*.log", mode: 'copy'

  input:
  tuple sampleID, file(cel_list_file) from cel_list

  output:
  tuple sampleID, file("lrr_baf.txt") into lrrbaf, dummy_lrrbaf

  script:
  """

  /apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-probeset-summarize --cdf-file ${params.snp6chip} --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch ${params.hapmap_norm_target} --out-dir ${sample_outputdir}_tmp  --cel-files ${cel_list_file}

  ${params.norm_geno_clust} ${params.genoclust} ${sample_outputdir}_tmp/quant-norm.pm-only.med-polish.expr.summary.txt -locfile ${params.gw6_pfb_file} -out lrr_baf.txt

  rm -rf ${sample_outputdir}_tmp/quant-norm.pm-only.med-polish.expr.summary.txt ${sample_outputdir}_tmp/quant-norm.pm-only.med-polish.expr.report.txt

  """

}


// Step 3: ASCAT
process hCNV_ASCAT {
  tag "sampleID"
  label 'short_mem'
  label 'ASCAT_and_annot'

  publishDir "${sample_outputdir}_tmp", pattern: "*.txt", mode: 'copy'
  publishDir "${sample_outputdir}_tmp", pattern: "*.png", mode: 'copy'
  publishDir "${sample_outputdir}_tmp", pattern: "*.Rout", mode: 'copy'
  publishDir "${sample_outputdir}_tmp", pattern: "*RData", mode: 'copy'

  input:
  tuple sampleID, file(lrr_baf) from lrrbaf
  tuple sampleID, file(gender) from gender1

  output:
  tuple sampleID, file("*segments_raw.txt") into raw_seg, dummy_raw_seg
  tuple sampleID, file("*ploidy.txt") into ploidy, dummy_ploidy
  file("*txt")
  file("*png")
  file("*Rout")
  file("*RData")

  shell:
  '''

  cut -f1-5 !{lrr_baf} > lrr_baf1.txt

  R CMD BATCH --slave "--args !{params.SNPpos} !{params.GC}" !{params.tumor_lrrbaf_ASCAT}

  '''

}


// Step 4: Annotation
process hCNV_annotation {
  tag "sampleID"
  label 'short_mem'
  label 'ASCAT_and_annot'

  publishDir "${sample_outputdir}_tmp", pattern: "*segments*.txt", mode: 'copy'

  input:
  tuple sampleID, file(raw_segs) from raw_seg
  tuple sampleID, file(sample_ploidy) from ploidy
  tuple sampleID, file(sample_gender) from gender2

  output:
  tuple sampleID, file("*segments_raw.extend.ensgene_cnvbreak.txt") into dummy_annot
  file("*seg*txt")

  shell:
  '''

  !{params.extend_raw_seg} !{raw_segs} !{sample_ploidy} !{params.chr_arm} !{sample_gender}

  !{params.ensembl_cnv} !{sampleID}.segments_raw.extend.txt !{params.exp_mart_genes}

  !{params.msp_get} !{sampleID}.segments_raw.extend.txt > msp.txt

  msp=(cat msp.txt)


  R CMD BATCH --slave "--args !{sampleID}.segments_raw.extend.txt $msp ./ " !{params.plot_seg}
  
  rm -rf tmp.txt

  '''

}

// Step 5: Renaming sample directory
process hCNV_transfer_files {
  tag "sampleID"
  label 'vshort_mem'

  input:
  tuple sampleID, file(cellist) from dummy_cel_list
  tuple sampleID, file(genderfile) from dummy_gender
  tuple sampleID, file(lrrbaffile) from dummy_lrrbaf
  tuple sampleID, file(rawsegmentsfile) from dummy_raw_seg
  tuple sampleID, file(ploidyfile) from dummy_ploidy
  tuple sampleID, file(cnvbreakfile) from dummy_annot

  script:
  log.info "-----Re-naming sample output directory for ${sampleID}-----"
  """

  mv "${sample_outputdir}_tmp"/* "${sample_outdir}"

  """
  }


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Summary ~ ~ ~ ~ ~
workflow.onComplete = {
    println "--------------------------"
    println "Pipeline Execution Summary"
    println "--------------------------"
    println "Completed at: ${workflow.complete}"
    println "Duration    : ${workflow.duration}"
    println "exit status : ${workflow.exitStatus}"
    println "Success   	 : ${workflow.success}"
    println "Error       : ${workflow.errorMessage}"
    println "Report      : ${workflow.errorReport}"
    }
