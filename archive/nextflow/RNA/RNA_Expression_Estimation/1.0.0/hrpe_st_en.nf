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
    params.ref_prefix           = null // absolute path to xenome transcriptome index files with index file name prefix included
    params.rsem_ref_prefix      = null // absolute path to rsem index files with index file name prefix included
    params.hsa_accession        = null // absolute path to human unplaced contigs file
    params.ref_fa               = null // absolute path to reference genome used for alignment
    params.ref_flat             = null // absolute path to reference in flat file format
    params.ribo_intervals       = null // absolute path to ribosomal intervals file
    params.probes               = null // absolute path to probes bed file
    params.ctp_genes            = null // absolute path to CTP genes bed file
    params.classifier_table     = null // absolute path to file containing list of genes used to determine if models should be classified as lymphoma
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
        --ref_prefix              absolute path to xenome transcriptome index files with index file name prefix included
        --rsem_ref_prefix         absolute path to rsem index files with index file name prefix included
        --hsa_accession           absolute path to human unplaced contigs file
        --ref_fa                  absolute path to reference genome used for alignment
        --ref_flat                absolute path to reference in flat file format
        --ribo_intervals          absolute path to ribosomal intervals file
        --probes                  absolute path to probes bed file
        --ctp_genes               absolute path to CTP genes bed file
        --classifier_table        absolute path to file containing list of genes used to determine if models should be classified as lymphoma


    Optional:
        --tmpdir                Path to hold temporary files for software execution.
        --preserve_work         if 'yes' then don't delete the ./work dir at end

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
if ( ! params.ref_prefix ) {
    exit 1, "Parameter ERROR: absolute path to xenome transcriptome index files with index file name prefix included must be specified."
}
if ( ! params.rsem_ref_prefix ) {
    exit 1, "Parameter ERROR: absolute path to rsem index files with index file name prefix included must be specified."
}
if ( ! params.hsa_accession ) {
    exit 1, "Parameter ERROR: absolute path to human unplaced contigs file must be specified."
}
if ( ! params.ref_fa ) {
    exit 1, "Parameter ERROR: absolute path to reference genome used for alignment must be specified."
}
if ( ! params.ref_flat ) {
    exit 1, "Parameter ERROR: absolute path to reference in flat file format must be specified."
}
if ( ! params.ribo_intervals ) {
    exit 1, "Parameter ERROR: absolute path to ribosomal intervals file must be specified."
}
if ( ! params.probes ) {
    exit 1, "Parameter ERROR: absolute path to probes bed file must be specified."
}
if ( ! params.ctp_genes ) {
    exit 1, "Parameter ERROR: absolute path to CTP genes bed file must be specified."
}
if ( ! params.classifier_table ) {
    exit 1, "Parameter ERROR: ($params.targets_exons_picard) absolute path to file containing list of genes used to determine if models should be classified as lymphoma must be specified."
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
println "Pipeline:       ${workflow.manifest.name}"
println "Description:    ${workflow.manifest.description}"
if(workflow.revision) {
  println "Pipeline Release: ${workflow.revision}"
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
println "Parameters......"
println ".  FASTQ inputs:                  ${params.fastqInputs}"
println ".  Prefix for xenome index files: ${params.ref_prefix}"
println ".  Prefix for rsem index files:   ${params.rsem_ref_prefix}"
println ".  Human unplaced contigs file:   ${params.hsa_accession}"
println ".  Reference genome:              ${params.ref_fa}"
println ".  Reference flat file:           ${params.ref_flat}"
println ".  Ribosomal intervals files:     ${params.ribo_intervals}"
println ".  Probes bed file:               ${params.probes}"
println ".  CTP genes bed file:            ${params.ctp_genes}"
println ".  Gene list for classification:  ${params.classifier_table}"
println ".  TMP dir:                       ${params.tmpdir}"
println "Run Start Time:  ${workflow.start}"


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
process hrpe_st_en_qual_stat {
  tag "sampleID"
  label 'long_mem'
  label 'python_2_7_3'

  publishDir ".", pattern: "*fastq.gz_stat", mode: 'copy'


  input:
  tuple sampleID, read1, read2 from sample_fastqs_ch

  output:
  tuple sampleID, file("${sampleID}_R{1,2}*filtered_trimmed") into trimmed_fastq
  file "*.fastq.gz_stat"
  tuple sampleID, file("*.fastq.gz_stat") into fq_stats,dummy_fq_stats

  script:
  log.info "-----Qual_Stat running on ${sampleID}-----"
  """

  python ${params.filter_trim} -M ${params.min_pct_hq_reads} ${read1} ${read2}

  """
  }


// Step 2a: Xenome
process hrpe_st_en_xenome_classification_a {
  tag "sampleID"
  label 'long_high_mem'
  label 'xenome'

  publishDir ".", pattern: "*.txt", mode: 'copy'

  input:
  tuple sampleID, file(trimmed) from trimmed_fastq

  output:
  tuple sampleID, file("human*{1,2}.fastq") into xenome_classified_fastq
  file "*.txt"
  tuple sampleID, file("*.txt") into xenome_stats, xenome_stats2, dummy_xenome_stats

  script:
  log.info "-----Xenome running on ${sampleID}-----"
  """

  /xenome-1.0.1-r/xenome classify -T 12 -P ${params.ref_prefix} --pairs --host-name mouse --graft-name human -i ${trimmed[0]} -i ${trimmed[1]} > ${sampleID}_xenome_stats.txt

  rm -rf *both*fastq* *mouse*fastq* *neither*fastq* *ambiguous*fastq*

  """
  }


// Step 2b: Sorting Xenome fastq files
process hrpe_st_en_xenome_classification_b {
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


// Step 3: RSEM
process hrpe_st_en_rsem_alignment_exp{
  tag "sampleID"
  label 'long_high_mem'
  label 'rsem'

  publishDir ".", pattern: "*stats", mode: 'copy'
  publishDir ".", pattern: "*results*", mode: 'copy'

  input:
  tuple sampleID, file(trimmed_sorted) from sorted_xenome_classified_fastq

  output:
  file "*stats"
  file "*results*"
  tuple sampleID, file("*aln.stats") into rsem_stats
  tuple sampleID, file("*genome.bam") into genome_sorted_bam
  tuple sampleID, file("*genes.results") into results_genes, dummy_results_genes
  tuple sampleID, file("*isoforms.results") into results_isoforms, dummy_results_isoforms

  script:
  log.info "-----Genome alignment running on ${sampleID}-----"

  if (params.reads == "stranded"){
    prob="--forward-prob 0"
  }
  if (params.reads == "non_stranded"){
    prob="--forward-prob 0.5"
  }

  """

  rsem-calculate-expression -p 8 --phred33-quals --seed-length ${params.seed_length} $prob --time --output-genome-bam ${params.aligner} --paired-end ${trimmed_sorted[0]} ${trimmed_sorted[1]} ${params.rsem_ref_prefix} ${sampleID} 2> ${sampleID}_rsem_aln.stats

  """
  }


// Step 4: Add gene name and normalization
process hrpe_st_en_gene_name_norm{
  tag "sampleID"
  label 'vshort_mem'
  label 'perl_R'

  publishDir ".", pattern: "*Normalized", mode: 'copy'


  input:
  tuple sampleID, file(genes), file(isoforms) from results_genes.join(results_isoforms)

  output:
  tuple sampleID, file("*genes.results.Normalized") into norm_gene_results, dummy_norm_gene_results
  file "*withGeneName"
  file "*.Normalized"

  script:
  log.info "-----Gene name addition and normalization running on ${sampleID}-----"

  """

  perl ${params.gene_name_norm} \
  -i1 ${genes} \
  -i2 ${isoforms} \
  -a1 ${params.hsa_accession}

  """
  }


Channel.of( sampleID, fqR1p, fqR2p )
    .toList()
    .set { sample_fastqs_ch2 }


// Step 5: Get Read Group Information
process hrpe_st_en_read_group {
  tag "sampleID"
  label 'vshort_mem'
  label 'python_2_7_3'

  input:
  tuple sampleID, file(read1), file(read2) from sample_fastqs_ch2

  output:
  tuple sampleID, file("*.txt") into read_grp1, read_grp2

  script:
  log.info "-----Read group information determination running on ${sampleID}-----"

  """

  python ${params.read_grp_det} -p -o ${sampleID}_read_group.txt ${read1}

  """
  }


// Step 6: Picard Alignment Metrics
process hrpe_st_en_picard_aln_metrics {
  tag "sampleID"
  label 'med_mem'
  label 'picard_metrics'

  publishDir ".", pattern: "*txt", mode: 'copy'
  publishDir ".", pattern: "*pdf", mode: 'copy'

  input:
  tuple sampleID, file(gen_sort_bam), file(read_grp) from genome_sorted_bam.join(read_grp1)

  output:
  tuple sampleID, file("*reorder_sort.bam") into reordered_sorted_bam
  file "*picard_aln_metrics.txt"
  tuple sampleID, file("*metrics.txt") into picard_metrics, dummy_picard_metrics
  file "*coverage_vs_transcript_plot.pdf"

  script:
  log.info "-----Picard alignment metrics running on ${sampleID}-----"

  if (params.reads == "stranded"){
      strand="STRAND=SECOND_READ_TRANSCRIPTION_STRAND"
    }
  if (params.reads == "non_stranded"){
      strand="STRAND=NONE"
    }

  """

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar AddOrReplaceReadGroups \
  INPUT=${gen_sort_bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_group.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_grp) \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar ReorderSam \
  INPUT=${sampleID}_genome_bam_with_read_group.bam \
  OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  REFERENCE=${params.ref_fa} \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar SortSam \
  SO=coordinate \
  INPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  OUTPUT=${sampleID}_reorder_sort.bam \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar CollectRnaSeqMetrics \
  I=${sampleID}_reorder_sort.bam \
  O=${sampleID}_picard_aln_metrics.txt \
  REF_FLAT=${params.ref_flat} \
  RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
  $strand \
  CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf

  """
  }


// Step 7a: GATK Coverage Stats
process hrpe_st_en_gatk_cov_stats_a {
  tag "sampleID"
  label 'long_mem'
  label 'gatk'


  input:
  tuple sampleID, file(sorted_bam) from reordered_sorted_bam


  output:
  tuple sampleID, file("*gatk_temp3*") into gatk_stat3
  tuple sampleID, file("*gatk_temp6*") into gatk_stat6

  script:
  log.info "-----GATK coverage stats, part 1 running on ${sampleID}-----"

  """

  samtools index ${sorted_bam}


  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp1.txt \
  -I ${sorted_bam} \
  -L  ${params.probes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS


  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp4.txt \
  -I ${sorted_bam} \
  -L ${params.ctp_genes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS


  ${params.gatk_form} ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt ${params.probes}

  ${params.gatk_form} ${sampleID}_gatk_temp4.txt ${sampleID}_gatk_temp5.txt ${sampleID}_gatk_temp6.txt ${params.ctp_genes}

  """
  }



// Step 7b: GATK Coverage Stats

process hrpe_st_en_gatk_stats_b {
  tag "sampleID"
  label 'long_mem'
  label 'python_2_7_3'

  publishDir ".", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(gatk_temp3), file(gatk_temp6) from gatk_stat3.join(gatk_stat6)


  output:
  file "*exome_interval_avg_median_coverage.bed"
  file "*CCP_interval_avg_median_coverage.bed"
  tuple sampleID, file("*CP_interval_avg_median_coverage.bed") into dummy_avg_med_cov_bed

  script:
  log.info "-----GATK coverage stats, part 2 running on ${sampleID}-----"

  """

  python ${params.cov_calc} ${sampleID}_gatk_temp3.txt ${sampleID}_exome_interval_avg_median_coverage.bed

  python ${params.cov_calc} ${sampleID}_gatk_temp6.txt ${sampleID}_CCP_interval_avg_median_coverage.bed

  """
  }


// Step 8: Classifier and Coverage
process hrpe_st_en_classifier_coverage{
  tag "sampleID"
  label 'short_mem'
  label 'class_cov'

  publishDir ".", pattern: "*classification", mode: 'copy'

  input:
  tuple sampleID, file(norm_genes) from norm_gene_results

  output:
  file "*classification"
  tuple sampleID, file("*classification") into classified


  script:
  log.info "-----Classifier and coverage running on ${sampleID}-----"
  """

  python ${params.lymph_class} -o ${sampleID}_classification ${norm_genes} ${params.classifier_table} ${sampleID}

  """
  }


// Step 9: Summary metrics
process hrpe_st_en_summary_metrics {
  tag "sampleID"
  label 'vshort_mem'
  label 'summ_metrics'

  publishDir ".", pattern: "*stats.txt", mode: 'copy'

  input:
  tuple sampleID, file(fastq_stat), file(xenome_stat), file(rsem_stat), file(picard_aln_met) from fq_stats.join(xenome_stats2).join(rsem_stats).join(picard_metrics)


  output:
  file "*summary_stats.txt"
  tuple sampleID, file("*stats.txt") into dummy_stats

  script:
  log.info "-----Summary metrics running on ${sampleID}-----"
  """

  perl ${params.sum_mets} ${fastq_stat} ${xenome_stat} ${rsem_stat} ${picard_aln_met} > ${sampleID}_summary_stats.txt

  """
  }



// Step 10: Finalization .. do some small file manipulation wrapup tasks, and
// remove 99% of the (large) ./work dir, unless preserve_work = "yes".
// We use errorStrategy 'ignore' because the "rm -rf ./work" seems to remove almost 
// everything in ./work ok, but then returns a status code of 1 likely due to the 
// conflict of one subdir in ./work still in use by this pipeline step.
process hrpe_st_en_finalization{
  tag "sampleID"
  label 'finalize'

  errorStrategy 'ignore'

  input:
  tuple sampleID, file(fastq_stat) from dummy_fq_stats
  tuple sampleID, file(xenome_stats) from dummy_xenome_stats
  tuple sampleID, file(gene_results) from dummy_results_genes
  tuple sampleID, file(isoform_results) from dummy_results_isoforms
  tuple sampleID, file(norm_gene_results) from dummy_norm_gene_results
  tuple sampleID, file(picard_metrics) from dummy_picard_metrics
  tuple sampleID, file(avg_med_coverage) from dummy_avg_med_cov_bed
  tuple sampleID, file(sum_stats) from dummy_stats

  script:
  log.info "-----Finalizing ${sampleID} ----- preserve_work = ${params.preserve_work}"
  """
  cd ${workflow.launchDir}
  mv  pipeline_running.txt  pipeline_completed.txt
  touch  pipeline_completed.txt
  mv *.fastq.gz_stat ${sampleID}_fastqs_stat.txt

  if [ "${params.preserve_work}" == "no" ]; then
    # this should be last ... raises an error
    rm  -rf  ./work
  fi

  """
  }
