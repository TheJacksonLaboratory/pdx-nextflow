#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show()



//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Configurable variable parameters specific to individual runs:
    params.fastqInputs          = null // absolute path to input fastq(s).
    params.tmpdir               = "$TMPDIR"
    params.refprefix            = null // absolute path to bwa index files with index file name prefix included
    params.ref_fa_bwa           = null // absolute path to bwa index .fa file
    params.ref_fa               = null // absolute path to reference genome to be used for alignment
    params.gold_std_indels      = null // absolute path to gold standard indels vcf file
    params.known_indels         = null // absolute path to known indels vcf file
    params.targets_gatk         = null // absolute path to CTP targets bed file for gatk
    params.dbsnp                = null // absolute path to dbSNP vcf file
    params.targets_picard       = null // absolute path to CTP targets bed file for picard
    params.targets_exons_picard = null // absolute path to CTP targets exons bed file for picard
    params.msisensor_model      = null // absolute path to model files required by MSIsensor2
    params.fa2bit               = null // absolute path to reference genome file in 2bit format
    params.ensembl_transcript   = null // absolute path to ensembl transcript list for snpEff
    params.dbNSFP               = null // absolute path to dbNSFP file
    params.Cosmic_newer         = null // absolute path to Cosmic variants file (v80)
    params.Cosmic_older         = null // absolute path to Cosmic variants file (v75)
    params.snpEff_data          = null // absolute path to snpEff data
    params.hgvs_data            = null // absolute path to hgvs data
    params.rec_var              = null // absolute path to CTP recurring variants file
    params.bins_ctpcoverage     = null // absolute path to CTP coverage bed file
    params.ctp_genes            = null // absolute path to CTP genes bed file
    params.preserve_work        = "no"

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
        --fastqInputs             absolute path to input fastq(s).
        --refprefix               absolute path to bwa index files with index file name prefix included
        --ref_fa_bwa              absolute path to bwa index .fa file
        --ref_fa                  absolute path to reference genome to be used for alignment
        --gold_std_indels         absolute path to gold standard indels
        --known_indels            absolute path to known indels
        --targets_gatk            absolute path to CTP targets bed file for gatk
        --dbsnp                   absolute path to dbSNP vcf file
        --targets_picard          absolute path to CTP targets bed file for picard
        --targets_exons_picard    absolute path to CTP targets exons bed file for picard
        --msisensor_model         absolute path to model files required by MSIsensor2
        --fa2bit                  absolute path to reference genome file in 2bit format
        --ensembl_transcript      absolute path to ensembl transcript list for snpEff
        --dbNSFP                  absolute path to dbNSFP file
        --Cosmic_newer            absolute path to newer (v80) Cosmic variants vcf file
        --Cosmic_older            absolute path to older (v75) Cosmic variants vcf file
        --snpEff_data             absolute path to snpEff data
        --hgvs_data               absolute path to hgvs data
        --rec_var                 absolute path to CTP recurring variants file
        --bins_ctpcoverage        absolute path to CTP coverage bed file
        --ctp_genes               absolute path to CTP genes bed file



    Optional:
        --tmpdir                Path to hold temporary files for software execution.

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
////// Required parameters \\\\\\
if ( ! params.fastqInputs ) {
    exit 1, "Parameter ERROR: fastqInputs ($params.fastqInputs) must be the absolute path to fastq file(s) with _R1 and _R2 fastq filenames."
}
if ( ! params.refprefix ) {
    exit 1, "Parameter ERROR: absolute path to bwa index files with index file name prefix included."
}
if ( ! params.ref_fa_bwa ) {
    exit 1, "Parameter ERROR: absolute path to absolute path to bwa index .fa file must be specified."
}
if ( ! params.ref_fa ) {
    exit 1, "Parameter ERROR: absolute path to reference genome to be used for alignment must be specified."
}
if ( ! params.gold_std_indels ) {
    exit 1, "Parameter ERROR: ($params.gold_std_indels) absolute path to gold standard indels file must be specified."
}
if ( ! params.known_indels ) {
    exit 1, "Parameter ERROR: ($params.known_indels) absolute path to known indels file must be specified."
}
if ( ! params.targets_gatk ) {
    exit 1, "Parameter ERROR: ($params.targets_gatk) absolute path to CTP targets bed file for gatk must be specified."
}
if ( ! params.dbsnp ) {
    exit 1, "Parameter ERROR: ($params.dbsnp) absolute path to dbSNP file must be specified."
}
if ( ! params.targets_picard ) {
    exit 1, "Parameter ERROR: ($params.targets_picard) absolute path to CTP targets bed file for picard must be specified."
}
if ( ! params.targets_exons_picard ) {
    exit 1, "Parameter ERROR: ($params.targets_exons_picard) absolute path to CTP targets exons bed file for picard must be specified."
}
if ( ! params.msisensor_model ) {
    exit 1, "Parameter ERROR: ($params.msisensor_model) absolute path to model files required by MSIsensor2 must be specified."
}
if ( ! params.fa2bit ) {
    exit 1, "Parameter ERROR: ($params.fa2bit) absolute path to reference genome file in 2bit format must be specified."
}
if ( ! params.ensembl_transcript ) {
    exit 1, "Parameter ERROR: ($params.ensembl_transcript) absolute path to ensembl transcript list for snpEff must be specified."
}
if ( ! params.dbNSFP ) {
    exit 1, "Parameter ERROR: ($params.dbNSFP) absolute path to dbNSFP file must be specified."
}
if ( ! params.Cosmic_newer ) {
    exit 1, "Parameter ERROR: ($params.Cosmic_newer) absolute path to newer (v80) Cosmic variants vcf file must be specified."
}
if ( ! params.Cosmic_older ) {
    exit 1, "Parameter ERROR: ($params.Cosmic_older) absolute path to older (v75) Cosmic variants vcf file must be specified."
}
if ( ! params.snpEff_data ) {
    exit 1, "Parameter ERROR: ($params.snpEff_data) absolute path to snpEff data must be specified."
}
if ( ! params.hgvs_data ) {
    exit 1, "Parameter ERROR: ($params.hgvs_data) absolute path to hgvs data must be specified."
}
if ( ! params.rec_var ) {
    exit 1, "Parameter ERROR: ($params.rec_var) absolute path to CTP recurring variants file must be specified."
}
if ( ! params.bins_ctpcoverage ) {
    exit 1, "Parameter ERROR: ($params.bins_ctpcoverage) absolute path to CTP coverage bed file must be specified."
}
if ( ! params.ctp_genes ) {
    exit 1, "Parameter ERROR: ($params.ctp_genes) absolute path to CTP genes bed file must be specified."
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
println "Pipeline:             ${workflow.manifest.name}"
println "Description:          ${workflow.manifest.description}"
if(workflow.revision) {
    println "Pipeline Release:     ${workflow.revision}"
}
println "Run Name:       ${workflow.runName}"
println "User:           ${workflow.userName}"
println "Config Profile: ${workflow.profile}"
println "Config Files:   ${workflow.configFiles}"
println "Command Line:   ${workflow.commandLine}"
println "Nextflow Info:  v. ${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
println "Launch dir:     ${workflow.launchDir}"
println "Working dir:    ${workflow.workDir}"
println "Workflow dir:   ${workflow.projectDir}"

// Pipeline Params:
println "Parameters......"
println ".  FASTQ inputs:                ${params.fastqInputs}"
println ".  Prefix for bwa index files:  ${params.refprefix}"
println ".  bwa index .fa file:          ${params.ref_fa_bwa}"
println ".  Reference genome used for alignment:   ${params.ref_fa}"
println ".  Gold Standard Indels:        ${params.gold_std_indels}"
println ".  Known Indels:                ${params.known_indels}"
println ".  CTP targets bed file for gatk:         ${params.targets_gatk}"
println ".  dbSNP vcf file:                        ${params.dbsnp}"
println ".  CTP targets bed file for picard:       ${params.targets_picard}"
println ".  CTP targets exons bed file for picard: ${params.targets_exons_picard}"
println ".  Model files required by MSIsensor2:    ${params.msisensor_model}"
println ".  Reference genome file in 2bit format:  ${params.fa2bit}"
println ".  ensembl transcript list for snpEff:    ${params.ensembl_transcript}"
println ".  dbNSFP file:                           ${params.dbNSFP}"
println ".  Cosmic vcf file (newer):     ${params.Cosmic_newer}"
println ".  Cosmic vcf file (older):     ${params.Cosmic_older}"
println ".  snpEff data:                 ${params.snpEff_data}"
println ".  HGVS data:                   ${params.hgvs_data}"
println ".  CTP recurring variants file: ${params.rec_var}"
println ".  CTP coverage bed file:       ${params.bins_ctpcoverage}"
println ".  CTP genes bed file:          ${params.ctp_genes}"
println ".  TMP dir:                     ${params.tmpdir}"
println "Run Start Time:                 ${workflow.start}"


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def tmpdir    = params.tmpdir


def fqPairRE = ~/_R1.*\..*f[ast]*q.*$/

fqin = params.fastqInputs.tokenize(",")
fqR1 = file(fqin[0])


// Need to specify absolute path to fastq file(s) because relative paths cause symlink breakage
fqR1p = fqR1.toAbsolutePath().toString()
fqR2p = file(fqin[1]).toAbsolutePath().toString()

def sampleID = ( fqR1.name - fqPairRE )


//~~~~~~~~~~ Initial Channel of SampleID and Fastq Data ~~~~
Channel.of( sampleID, fqR1p, fqR2p )
       .toList()
       .set { sample_fastqs_ch }


// ~~~~~~~~~~ publishDir will be current working dir (.) ~~~~~


// Step 1: Qual_Stat
process hctp_qual_stat {
  tag "sampleID"
  label 'long_mem'
  label 'python_2_7_3'

  publishDir ".", pattern: "*fastqs_stat.txt", mode: 'copy'


  input:
  tuple sampleID, read1, read2 from sample_fastqs_ch

  output:
  tuple sampleID, file("${sampleID}_R{1,2}*filtered_trimmed") into trimmed_fastq
  file "*fastqs_stat.txt"
  tuple sampleID, file("*fastqs_stat.txt") into fq_stats,dummy_fq_stats


  script:
  log.info "-----Qual_Stat running on ${sampleID}-----"
  """

  python ${params.filter_trim} -M ${params.min_pct_hq_reads} ${read1} ${read2}

  mv *.fastq.gz_stat ${sampleID}_fastqs_stat.txt

  """
  }


// Step 2a: Xenome
process hctp_xenome_classification_a {
  tag "sampleID"
  label 'long_high_mem'
  label 'xenome'

  publishDir ".", pattern: "*.txt", mode: 'copy'

  input:
  tuple sampleID, file(trimmed) from trimmed_fastq

  output:
  tuple sampleID, file("human*{1,2}.fastq") into xenome_classified_fastq
  file "*.txt"
  tuple sampleID, file("*.txt") into xenome_stats,dummy_xenome_stats

  script:
  log.info "-----Xenome running on ${sampleID}-----"
  """

  /xenome-1.0.1-r/xenome classify -T 12 -P ${params.refprefix} --pairs --host-name mouse --graft-name human -i ${trimmed[0]} -i ${trimmed[1]} > ${sampleID}_xenome_stats.txt


  rm -rf *both*fastq* *mouse*fastq* *neither*fastq* *ambiguous*fastq*

  """
  }


// Step 2b: Sorting Xenome fastq files
process hctp_xenome_classification_b {
  tag "sampleID"
  label 'long_high_mem'
  label 'fq_tools'

  input:
  tuple sampleID, file(trimmed_hsa) from xenome_classified_fastq
  tuple sampleID, file(xen_stats) from xenome_stats

  output:
  tuple sampleID, file("*sorted_human_{1,2}.fastq") into sorted_xenome_classified_fastq

  script:
  log.info "-----Sorting human reads from Xenome step running on ${sampleID}-----"
  """

  fastq-sort --id ${trimmed_hsa[0]} > ${sampleID}_sorted_human_1.fastq
  fastq-sort --id ${trimmed_hsa[1]} > ${sampleID}_sorted_human_2.fastq

  """
  }


// Step 3: Get Read Group Information and BWA-MEM Alignment
process hctp_bwa_mem {
  tag "sampleID"
  label 'long_high_mem'
  label 'bwakit'

  input:
  tuple sampleID, file(trimmed_sorted) from sorted_xenome_classified_fastq


  output:
  tuple sampleID, file("*aln.bam") into bwa_mem_out,dummy_bwa_mem_out


  script:
  log.info "-----BWA-MEM Alignment running on ${sampleID}-----"
  """

  python ${params.read_grp} ${trimmed_sorted[0]} -o ${sampleID}_read_group.txt

  rg=\$(cat ${sampleID}_read_group.txt)

  run-bwamem -t ${task.cpus} -R \${rg} -o ${sampleID} -H ${params.ref_fa_bwa} ${trimmed_sorted[0]} ${trimmed_sorted[1]} | sh

  rm -rf *read_group.txt *merged.fastq *temp.sam *out

  """
  }


// Step 4: Variant Pre-Processing 1
process hctp_variant_preproc_1 {
  tag "sampleID"
  label 'short_mem'
  label 'picard_2_8_1'

  publishDir ".", pattern: "*.dat", mode: 'copy'
  publishDir ".", pattern: "*.bam*", mode: 'copy'

  input:
  tuple sampleID, file(init_bam) from bwa_mem_out


  output:
  tuple sampleID, file("*dedup.bam") into bam_dedup1
  tuple sampleID, file("*dedup.bai") into bam_dedup2
  tuple sampleID, file("*_dup_metrics.dat") into picard_metrics,dummy_picard_metrics
  file "*_dup_metrics.dat"
  file "*bam"

  script:
  log.info "-----Variant pre-processing 1 running on ${sampleID}-----"
  """

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar SortSam \
  SO=coordinate \
  INPUT=${init_bam} \
  OUTPUT=${sampleID}_sorted.bam \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar MarkDuplicates \
  I=${sampleID}_sorted.bam \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.dat \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT

  """
  }


// Step 5: Variant Pre-Processing 2
process hctp_variant_preproc_2 {
  tag "sampleID"
  label 'long_mem'
  label 'gatk3'

  publishDir ".", pattern: "*realigned_BQSR.bam", mode: 'copy'

  input:
  tuple sampleID, file(dedup_bam) from bam_dedup1
  tuple sampleID, file(dedup_bai) from bam_dedup2

  output:
  tuple sampleID, file("*realigned_BQSR.bam") into realigned_bam1, realigned_bam2, realigned_bam3, realigned_bam4, realigned_bam5, realigned_bam6
  tuple sampleID, file("*bai") into realigned_bam_index1, realigned_bam_index2, realigned_bam_index3
  file "*realigned_BQSR.bam"

  script:
  log.info "-----Variant pre-processing 2 running on ${sampleID}-----"
  """

  java -Djava.io.tmpdir=$TMPDIR -Xmx2g -jar /GenomeAnalysisTK.jar \
  -I ${dedup_bam} \
  -R ${params.ref_fa} \
  -T RealignerTargetCreator \
  -o aligner.intervals \
  -known ${params.gold_std_indels} \
  -known ${params.known_indels} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -L ${params.targets_gatk}

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /GenomeAnalysisTK.jar \
  -I ${dedup_bam} \
  -R ${params.ref_fa} \
  -T IndelRealigner \
  -targetIntervals aligner.intervals \
  -o realigned.bam \
  -known ${params.gold_std_indels} \
  -known ${params.known_indels} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -L ${params.targets_gatk}

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -I realigned.bam \
  -R ${params.ref_fa} \
  -knownSites ${params.dbsnp} \
  -knownSites ${params.gold_std_indels} \
  -knownSites ${params.known_indels} \
  -o recal_data.grp \
  --disable_auto_index_creation_and_locking_when_reading_rods

  java -Xmx4g -jar /GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${params.ref_fa} \
  -I realigned.bam \
  -BQSR recal_data.grp \
  -o ${sampleID}_realigned_BQSR.bam \
  --disable_auto_index_creation_and_locking_when_reading_rods

  samtools index ${sampleID}_realigned_BQSR.bam

  rm -rf aligner.intervals realigned.bam realigned.bai recal_data.grp

  """
  }


// Step 6: Variant Pre-Processing 3
process hctp_variant_preproc_3 {
  tag "sampleID"
  label 'vshort_mem'
  label 'picard_1_95_python_2_7_3'

  publishDir ".", pattern: "*.txt", mode: 'copy'

  input:
  tuple sampleID, file(bam_realigned) from realigned_bam1

  output:
  tuple sampleID, file("*_CoverageMetrics.txt") into cov_met,dummy_cov_met
  file "*_CoverageMetrics.txt"

  script:
  log.info "-----Variant pre-processing 3 running on ${sampleID}-----"
  """

  sampleID=\$(basename $bam_realigned _realigned_BQSR.bam)

  java -Djava.io.tmpdir=$TMPDIR -jar -Xmx2g /picard-tools-1.95/CalculateHsMetrics.jar \
  TARGET_INTERVALS= ${params.targets_picard} \
  BAIT_INTERVALS= ${params.targets_exons_picard} \
  REFERENCE_SEQUENCE=${params.ref_fa} \
  INPUT=${bam_realigned} \
  OUTPUT=${sampleID}_CoverageMetrics.txt \
  VALIDATION_STRINGENCY=LENIENT

  python ${params.filt_dna_cov} ${sampleID}_CoverageMetrics.txt

  """
  }


// Step 7: MSIsensor2
process hctp_msisensor2 {
  tag "sampleID"
  label 'long_mem'
  label 'gcc_zlib_glibc_msisensor2'

  publishDir ".", pattern: "*sensor*", mode: 'copy'

  input:
  tuple sampleID, file(bam_realigned) from realigned_bam2
  tuple sampleID, file(bam_realigned_index) from realigned_bam_index1

  output:
  tuple sampleID, file("*msisensor") into msisensor2out,dummy_msisensor2out
  file "*msi*"

  script:
  log.info "-----MSIsensor2 running on ${sampleID}-----"
  """

  mkdir models

  cp -r ${params.msisensor_model} models

  msisensor2 msi -M models/models_hg38 -t ${bam_realigned} -o ${sampleID}_msisensor

  """
  }


// Step 8: Variant Calling Mutect2
process hctp_var_call_mutect2 {
  tag "sampleID"
  label 'vlong_mem'
  label 'mutect2'

  publishDir ".", pattern: "*.vcf*", mode: 'copy'

  input:
  tuple sampleID, file(bam_realigned) from realigned_bam3

  output:
  tuple sampleID, file("*final_mutect_snp_indel_filtered.vcf") into mutect_snp_indel_filtered, dummy_snp_indel_filt
  file "*final*filtered.vcf*"
  file "*intermed.vcf*"
  tuple sampleID, file("*_intermed.vcf.gz") into dummy_mut_snp_indel

  script:
  log.info "-----Variant Calling Mutect2 running on ${sampleID}-----"
  """


  java -Xmx4g -jar /gatk/build/libs/gatk-package-4.0.5.1-local.jar \
  GetSampleName \
  -I ${bam_realigned} \
  -O tumor_SN.txt

  tumorName=\$(cat tumor_SN.txt)

  java -Djava.io.tmpdir=$TMPDIR -Xmx120g -jar /gatk/build/libs/gatk-package-4.0.5.1-local.jar \
  Mutect2 \
  -R ${params.ref_fa} \
  -I ${bam_realigned} \
  -tumor \${tumorName} \
  -O ${sampleID}_intermed.vcf \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --dont-use-soft-clipped-bases \
  --genotype-germline-sites \
  --sample-ploidy 4 \
  --annotation QualByDepth \
  --annotation RMSMappingQuality \
  --annotation FisherStrand \
  --annotation MappingQualityRankSumTest \
  --annotation ReadPosRankSumTest \
  -L ${params.targets_gatk}

  bgzip ${sampleID}_intermed.vcf

  tabix ${sampleID}_intermed.vcf.gz

  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /gatk/build/libs/gatk-package-4.0.5.1-local.jar \
  FilterMutectCalls \
  --variant ${sampleID}_intermed.vcf.gz \
  --output ${sampleID}_final_mutect_snp_indel_filtered.vcf \
  --min-base-quality-score 20 \
  --dont-use-soft-clipped-bases true \
  --unique-alt-read-count 5 \
  -stand-call-conf 30 \
  -L ${params.targets_gatk}

  rm -rf tumor_SN.txt

  """
  }


//Step 9: Variant Filtration
process hctp_var_filt{
  tag "sampleID"
  label 'vshort_mem'
  label 'variant_filtration'

  publishDir ".", pattern: "*filt*.vcf", mode: 'copy'

  input:
  tuple sampleID, file(mutect2_filtered) from mutect_snp_indel_filtered

  output:
  tuple sampleID, file("*variants.DPfiltered.vcf") into DP_filtered_var1, DP_filtered_var2, dummy_DP_filt
  file "*filtered.vcf"


  script:
  log.info "-----Variant Filtration running on ${sampleID}-----"
  """

  python ${params.AD_min_AF} ${mutect2_filtered} ${sampleID}_mutect_snp_indel_filtered.vcf.DPfiltered.tmp.vcf 140

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar \
  annotate \
  -id ${params.dbsnp} ${sampleID}_mutect_snp_indel_filtered.vcf.DPfiltered.tmp.vcf > ${sampleID}_mutect_snp_indel_filtered.vcf.additionalfilters.tmp.vcf


  ${params.add_caller_gatk} \
  ${sampleID}_mutect_snp_indel_filtered.vcf.additionalfilters.tmp.vcf ${sampleID}_variants.DPfiltered.vcf

  rm -rf ${sampleID}_mutect_snp_indel_filtered.vcf.DPfiltered.tmp.vcf ${sampleID}_mutect_snp_indel_filtered.vcf.additionalfilters.tmp.vcf

  """
  }


// Step 10: Fix adjacent SNPs
process hctp_fix_adj_snps{
  tag "sampleID"
  label 'short_mem'
  label 'adj_snps_fix'

  publishDir ".", pattern: "*AdjSNP*", mode: 'copy'

  input:
  tuple sampleID, file(BQSR_bam_realigned), file(filt_var_DP) from realigned_bam4.join(DP_filtered_var1)

  output:
  tuple sampleID, file("*variant_fixAdjSNP.vcf") into fixed_adj_snps,dummy_fixed_adj_snps
  file "*.vcf"
  file "*.log.txt"

  script:
  log.info "-----Fix Adjacent SNPs running on ${sampleID}-----"
  """

  python ${params.adj_snps} \
  -v ${filt_var_DP} \
  -o ${sampleID}_variant_fixAdjSNP.vcf \
  1 ${BQSR_bam_realigned} ${params.fa2bit} 2>${sampleID}_AdjSNP.log.txt

  awk 'NF' ${sampleID}_variant_fixAdjSNP.vcf > temp

  mv temp ${sampleID}_variant_fixAdjSNP.vcf

  bgzip ${sampleID}_variant_fixAdjSNP.vcf

  tabix ${sampleID}_variant_fixAdjSNP.vcf.gz

  bcftools annotate \
  --output ${sampleID}_variant_fixAdjSNP.noIds.vcf.gz \
  --output-type z \
  --remove ID ${sampleID}_variant_fixAdjSNP.vcf.gz


  java -jar /snpEff_v4_3/snpEff/SnpSift.jar annotate \
  -id ${params.dbsnp} ${sampleID}_variant_fixAdjSNP.noIds.vcf.gz > ${sampleID}_variant_fixAdjSNP.vcf

  """
  }


// Step 11a: microIndel calling
process hctp_microindel_calling_a{
  tag "sampleID"
  label 'long_mem'
  label 'microIndel_call_a'

  publishDir ".", pattern: "*microIndels*vcf", mode: 'copy'

  input:
  tuple sampleID, file(BQSR_realigned_bam) from realigned_bam5
  tuple sampleID, file(bam_realigned_index) from realigned_bam_index2

  output:
  tuple sampleID, file("*microIndels_all.vcf") into all_microIndels
  file "*vcf"

  script:
  log.info "-----Microindel calling, part 1 running on ${sampleID}-----"
  """

  echo -e "${BQSR_realigned_bam}\t350\t${sampleID}" > ${sampleID}_pindel_config.txt

  pindel \
  --fasta ${params.ref_fa} \
  --config-file ${sampleID}_pindel_config.txt \
  -o ${sampleID}

  cat  ${sampleID}_D  ${sampleID}_SI >  ${sampleID}_DSI

  /pindel-0.2.5/pindel2vcf \
  -p  ${sampleID}_DSI \
  -r ${params.ref_fa} \
  -R hg38 \
  -d 20150925 \
  --max_size 50 \
  --vcf ${sampleID}_microIndels_all.vcf \
  -G \
  --het_cutoff 0.05

  """
  }

// Step 11b: microIndel calling
process hctp_microindel_calling_b {
  tag "sampleID"
  label 'long_mem'
  label 'microIndel_call_b'

  publishDir ".", pattern: "*.vcf", mode: 'copy'

  input:
  tuple sampleID, file(microIndels_all) from all_microIndels

  output:
  file "*.vcf"
  tuple sampleID, file("*microIndels.DPfiltered1.vcf") into DP_filt_indels,dummy_DP_filt_indels

  script:
  log.info "-----microIndel calling, part 2 running on ${sampleID}-----"
  """

  bedtools intersect \
  -header \
  -a ${microIndels_all} \
  -b  ${params.targets_gatk} \
  -f 1.0 > ${sampleID}_microIndels.raw.vcf

  python ${params.min_DP_filt} ${sampleID}_microIndels.raw.vcf ${sampleID}_microIndels.DPfiltered1.vcf

  touch  ${sampleID}_BP  ${sampleID}_D  ${sampleID}_DSI  ${sampleID}_INT  ${sampleID}_INT_final  ${sampleID}_INV  ${sampleID}_LI  ${sampleID}_RP  ${sampleID}_SI  ${sampleID}_TD

  rm -rf ${sampleID}_BP ${sampleID}_CloseEndMapped ${sampleID}_D ${sampleID}_DSI ${sampleID}_INT ${sampleID}_INT_final ${sampleID}_INV ${sampleID}_LI ${sampleID}_RP ${sampleID}_SI ${sampleID}_TD

  """
  }


// Step 12: microIndel filtration
process hctp_microindel_filtering{
  tag "sampleID"
  label 'vshort_mem'
  label 'microIndel_filt'

  publishDir ".", pattern: "*DPfiltered.vcf", mode: 'copy'

  input:
  tuple sampleID, file(indels_dp_filt) from DP_filt_indels

  output:
  tuple sampleID, file("*microIndels.DPfiltered.vcf") into DP_indels_filt1, DP_indels_filt2,dummy_DP_indels_filt
  file "*DPfiltered.vcf"

  script:
  log.info "-----microIndel filtration running on ${sampleID}-----"
  """

  python ${params.AD_min_AF} ${indels_dp_filt} ${sampleID}_microIndels.DPfiltered1.vcf.DPfiltered_microindels.tmp.vcf 140


  java -jar /snpEff_v4_3/snpEff/SnpSift.jar \
  annotate -id ${params.dbsnp} ${sampleID}_microIndels.DPfiltered1.vcf.DPfiltered_microindels.tmp.vcf > ${sampleID}_microIndels.DPfiltered1.vcf.additionalfilters.tmp.vcf

  ${params.add_caller_pindel} \
  ${sampleID}_microIndels.DPfiltered1.vcf.additionalfilters.tmp.vcf ${sampleID}_microIndels.DPfiltered.vcf

  rm -rf ${sampleID}_microIndels.DPfiltered1.vcf.DPfiltered_microindels.tmp.vcf ${sampleID}_microIndels.DPfiltered1.vcf.additionalfilters.tmp.vcf

  """
  }


// Step 13: Variant Annotation
process hctp_variant_annot{
  tag "sampleID"
  label 'short_mem'
  label 'variant_annot'

  publishDir ".", pattern: "*Annotated.*" , mode: 'copy'

  input:
  tuple sampleID, file (adj_snps_fixed), file(filtered_indels_DP) from fixed_adj_snps.join(DP_indels_filt1)

  output:
  file "*variants*Annotated.t*"
  tuple sampleID, file("*microIndels.Hardfiltered.Annotated.txt") into dummy_Hardfiltered_Annot

  shell:
  log.info "-----Variant annotation running on ${sampleID}-----"
  '''

  cat !{adj_snps_fixed} !{filtered_indels_DP} > !{sampleID}_temp.vcf

  vcf-sort !{sampleID}_temp.vcf > !{sampleID}_all_genes_variants_microindels.vcf


  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' !{sampleID}_all_genes_variants_microindels.vcf | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "PASS" --rmFilter "clustered_events" "( FILTER == 'clustered_events' )" | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;minDP' )" | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk' )" | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --rmFilter "clustered_events" "( FILTER == 'clustered_events;germline_risk;minDP' )" > !{sampleID}_all_genes_variants_microindels_filtered.vcf


  java -jar /snpEff_v4_3/snpEff/snpEff.jar eff -c !{params.snpEff_config} -v -lof -dataDir !{params.hgvs_data} -onlyTr !{params.ensembl_transcript} -hgvs GRCh38.84 -noStats !{sampleID}_all_genes_variants_microindels_filtered.vcf > !{sampleID}_all_genes_variants_microindels_snpEff.vcf


  java -jar /snpEff_v4_3/snpEff/SnpSift.jar dbnsfp -v -db !{params.dbNSFP} -noDownload -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AC,ExAC_AF !{sampleID}_all_genes_variants_microindels_snpEff.vcf > !{sampleID}_all_genes_variants_microindels_snpEff_snpSift.vcf


  cat !{sampleID}_all_genes_variants_microindels_snpEff_snpSift.vcf | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_all_genes_variants_microindels_snpEff_snpSift_onePerLine.vcf


  java -jar /snpEff_v4_3/snpEff/SnpSift.jar annotate -id !{params.Cosmic_newer} !{sampleID}_all_genes_variants_microindels_snpEff_snpSift_onePerLine.vcf > !{sampleID}_all_genes_variants_microindels_cosmicannotation.vcf


  java -jar /snpEff_v4_3/snpEff/SnpSift.jar extractFields !{sampleID}_all_genes_variants_microindels_cosmicannotation.vcf CHROM POS REF ALT ID FILTER DP_HQ ALT_AF "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_ExAC_AC" "dbNSFP_ExAC_AF" "CALLER" > !{sampleID}_variants_microIndels.DPfiltered.Annotated.tab


  python !{params.intergen_reg_names} -f !{sampleID}_variants_microIndels.DPfiltered.Annotated.tab -d


  cat !{sampleID}_variants_microIndels.DPfiltered.Annotated.tab | awk -F '\t' 'BEGIN {OFS="\t"} $6 == "FILTER" || $6 == "PASS" || $6 == "" || $6 == "."' > !{sampleID}_variants_microIndels.Hardfiltered.Annotated.txt

  rm -rf !{sampleID}_all_genes_variants_microindels.vcf !{sampleID}_all_genes_variants_microindels_filtered.vcf !{sampleID}_all_genes_variants_microindels_snpEff.vcf !{sampleID}_all_genes_variants_microindels_snpEff_snpSift.vcf !{sampleID}_all_genes_variants_microindels_snpEff_snpSift_onePerLine.vcf !{sampleID}_all_genes_variants_microindels_cosmicannotation.vcf !{sampleID}_all_genes_variants_microindels_cosmicannotation_beforeintegenic.vcf

  '''
  }


// Step 14: Tumor Mutation Burden Score
// High TMB cutoff is >=22
// High, moderate; FILTER=PASS; not in dbSNP

process hctp_tmb_score {
  tag "sampleID"
  label 'vshort_mem'
  label 'tmb_score'

  publishDir ".", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file (UG_vcf), file (pindel_vcf) from DP_filtered_var2.join(DP_indels_filt2)

  output:
  file "*.tab"
  file "*score"
  tuple sampleID, file("*HM.tab") into HM_TMB, dummy_HM_TMB
  tuple sampleID, file("*score") into dummy_TMB_score

  shell:
  log.info "-----TMB Score Calculation running on ${sampleID}-----"
  '''

  cat !{UG_vcf} !{pindel_vcf} > !{sampleID}_UG_pindel.vcf

  bedtools intersect -header -v -a !{sampleID}_UG_pindel.vcf -b !{params.rec_var}/*bed | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowAF" --rmFilter "PASS" 'ALT_AF[ANY] < 5' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "strandBias" --rmFilter "PASS" 'FS > 60' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowMQ" --rmFilter "PASS" 'MQ < 40' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowMQRankSum" --rmFilter "PASS" 'MQRankSum < -12.5' | \
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "lowReadPosRankSum" --rmFilter "PASS" 'ReadPosRankSum < -8' > !{sampleID}_nofalsePositives.tmp2.vcf

  java -jar /snpEff_v4_3/snpEff/snpEff.jar eff -v -dataDir !{params.snpEff_data} -lof -canon -hgvs hg38 -noStats !{sampleID}_nofalsePositives.tmp2.vcf > !{sampleID}_all_genes_variants_snpEff.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar dbnsfp -v -db !{params.dbNSFP} -noDownload -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AC,ExAC_AF !{sampleID}_all_genes_variants_snpEff.vcf > !{sampleID}_all_genes_variants_snpeff_snpsift.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar annotate -id !{params.Cosmic_older} !{sampleID}_all_genes_variants_snpeff_snpsift.vcf > !{sampleID}_all_genes_variants_cosmicannotation.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar filter --addFilter "PutativeGermline" --rmFilter "PASS" '(ALT_AF[ANY] >= 90 & dbNSFP_1000Gp3_AF[ANY] >= 0.0095) | ((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_1000Gp3_AF[ANY] >= 0.0095)| (ALT_AF[ANY] >= 90 & dbNSFP_ESP6500_AA_AF[ANY] >= 0.0095)| ((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ESP6500_AA_AF[ANY] >= 0.0095)|(ALT_AF[ANY] >= 90 & dbNSFP_ESP6500_EA_AF[ANY] >= 0.0095)|((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ESP6500_EA_AF[ANY] >= 0.0095)|(ALT_AF[ANY] >= 90 & dbNSFP_ExAC_AF[ANY] >= 0.0095)|((ALT_AF[ANY] >= 40 & ALT_AF[ANY] <= 60) & dbNSFP_ExAC_AF[ANY] >= 0.0095)' !{sampleID}_all_genes_variants_cosmicannotation.vcf > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag.vcf

  cat !{sampleID}_all_genes_variants_cosmicannotation_germlineflag.vcf | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline.vcf

  java -jar /snpEff_v4_3/snpEff/SnpSift.jar extractFields !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline.vcf CHROM POS REF ALT ID FILTER DP_HQ ALT_AF "LOF[*].NUMTR" "LOF[*].PERC" "EFF[*].GENE" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].RANK" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].CODING" "EFF[*].TRID" "dbNSFP_SIFT_score" "dbNSFP_SIFT_pred" "dbNSFP_Polyphen2_HDIV_score" "dbNSFP_MutationAssessor_score" "dbNSFP_phyloP100way_vertebrate" "dbNSFP_1000Gp3_AF" "dbNSFP_1000Gp3_AFR_AF" "dbNSFP_1000Gp3_EUR_AF" "dbNSFP_1000Gp3_AMR_AF" "dbNSFP_1000Gp3_EAS_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_ExAC_AC" "dbNSFP_ExAC_AF" "CALLER" > !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline_final.vcf


  header="chr\\tstart\\tend\\tCTPlength\\tHM"

  cat !{sampleID}_all_genes_variants_cosmicannotation_germlineflag_oneperline_final.vcf | grep "HIGH\\|MODERATE" | awk -F '\\t' '{ if($5 !~ "rs" && ($6==""||$6=="PASS"||$6==".")) print $1,$2-1,$2 }'| sort| uniq| tr ' ' '\\t' > !{sampleID}_count2 ; echo -e ${header} > !{sampleID}_HM.tab


  bedtools coverage -a !{params.bins_ctpcoverage} -b !{sampleID}_count2 | cut -f 1-5 >> !{sampleID}_HM.tab


  Rscript !{params.tmb_prog} !{sampleID}_HM.tab !{sampleID}_TMB.score

  '''
  }


// Step 15: Aggregate Stats
process hctp_agg_stats {
  tag "sampleID"
  label 'vshort_mem'
  label 'python_2_7_3'

  publishDir ".", pattern: "*stats.txt", mode: 'copy'

  input:
  tuple sampleID, file(filter_stat), file(duplicate_metrics), file(cov_metrics) from fq_stats.join(picard_metrics).join(cov_met)

  output:
  file "*_stats.txt"
  tuple sampleID, file("*summary_stats.txt") into dummy_aggregated_stats


  script:
  log.info "-----Aggregate Stats running on ${sampleID}-----"
  """

  python ${params.agg_stat} ${sampleID}_summary_stats.txt ${filter_stat} ${duplicate_metrics} ${cov_metrics}

  """
  }


// Step 16a: GATK Coverage Stats
process hctp_gatk_stats_a {
  tag "sampleID"
  label 'short_mem'
  label 'gatk_coverage_stats'

  input:
  tuple sampleID, file(BQSR_realigned_bam) from realigned_bam6
  tuple sampleID, file(bam_realigned_index) from realigned_bam_index3

  output:
  tuple sampleID, file("*temp3.txt") into gatk_stat

  script:
  log.info "-----GATK Coverage Stats, part 1 running on ${sampleID}-----"
  """

  java -Djava.io.tmpdir=$TMPDIR -Xmx12g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp1.txt \
  -I ${BQSR_realigned_bam} \
  -L  ${params.ctp_genes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable

  ${params.gatk_form} ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt  ${params.ctp_genes}

  """
  }


// Step 16b: GATK Coverage Stats
process hctp_gatk_stats_b {
  tag "sampleID"
  label 'short_mem'
  label 'python_2_7'

  publishDir ".", pattern: "*.bed", mode: 'copy'

  input:
  tuple sampleID, file(gatk_temp3) from gatk_stat

  output:
  tuple sampleID, file("*targetinterval_avg_median_coverage.bed") into dummy_avg_med_cov
  file "*.bed"

  script:
  log.info "-----GATK Coverage Stats, part 2 running on ${sampleID}-----"
  """

  python ${params.cov_calc} ${gatk_temp3} ${sampleID}_targetinterval_avg_median_coverage.bed

  """
  }



// Step 17: Finalization .. do some small file manipulation wrapup tasks, and
// remove 99% of the (large) ./work dir, unless preserve_work = "yes".
// We use errorStrategy 'ignore' because the "rm -rf ./work" seems to remove almost 
// everything in ./work ok, but then returns a status code of 1 likely due to the 
// conflict of one subdir in ./work still in use by this pipeline step.

// Step 17: Renaming sample directory
process hctp_finalization{
  tag "sampleID"
  label 'finalize'

  errorStrategy 'ignore'

  input:
  tuple sampleID, file(dummy_target_interval_avg_med_cov) from dummy_avg_med_cov
  tuple sampleID, file(dummy_Hardfilt_Annot) from dummy_Hardfiltered_Annot
  tuple sampleID, file(dummy_TMB_tab) from dummy_HM_TMB
  tuple sampleID, file(dummy_TMB_score) from dummy_TMB_score
  tuple sampleID, file(dummy_snp_indel_filt) from dummy_snp_indel_filt
  tuple sampleID, file(dummy_mut_snp_indel) from dummy_mut_snp_indel
  tuple sampleID, file(dummy_var_DPfilt) from dummy_DP_filt
  tuple sampleID, file(dummy_sum_stats) from dummy_aggregated_stats
  tuple sampleID, file(dummy_fq_stat) from dummy_fq_stats
  tuple sampleID, file(dummy_xenome_stat) from dummy_xenome_stats
  tuple sampleID, file(dummy_aln_bam) from dummy_bwa_mem_out
  tuple sampleID, file(dummy_picard_mets) from dummy_picard_metrics
  tuple sampleID, file(dummy_Coverage_Metrics) from dummy_cov_met
  tuple sampleID, file(dummy_msisensor) from dummy_msisensor2out
  tuple sampleID, file(dummy_variant_fixAdjSNP) from dummy_fixed_adj_snps
  tuple sampleID, file(dummy_microIndels_DPfiltered1) from dummy_DP_filt_indels
  tuple sampleID, file(dummy_microIndels_DPfiltered) from dummy_DP_indels_filt


  script:
  log.info "-----Finalizing ${sampleID} ----- preserve_work = ${params.preserve_work}"
  """
  cd ${workflow.launchDir}
  mv  pipeline_running.txt  pipeline_completed.txt
  touch  pipeline_completed.txt

  if [ "${params.preserve_work}" == "no" ]; then
    rm -f *microIndels_all.vcf
    rm -f *microIndels.DPfiltered1.vcf
    rm -f *microIndels.raw.vcf
    rm -f *.bam

    # this statement should be last... raises an error because ./work is in use by this step
    rm  -rf  ./work
  fi


  """
  }

