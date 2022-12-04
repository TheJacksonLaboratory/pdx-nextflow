#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from '../bin/shared/getLibraryId.nf'
include {param_log} from "${projectDir}/bin/log/wes"
include {RUN_START} from "${projectDir}/bin/shared/run_start"
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {QUALITY_STATISTICS} from '../modules/utility_modules/quality_stats'
include {XENOME_CLASSIFY} from   '../modules/xenome/xenome'
include {FASTQ_SORT as XENOME_SORT} from   '../modules/fastq-tools/fastq-tools_sort'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {BWA_MEM} from '../modules/bwa/bwa_mem'
include {PICARD_SORTSAM} from '../modules/picard/picard_sortsam'
include {PICARD_MARKDUPLICATES} from '../modules/picard/picard_markduplicates'
include {GATK_REALIGNERTARGETCREATOR} from '../modules/gatk/gatk_realignertargetcreator'
include {GATK_INDELREALIGNER} from '../modules/gatk/gatk_indelrealigner'
include {GATK_BASERECALIBRATOR} from '../modules/gatk/gatk_baserecalibrator'
include {GATK_PRINTREADS} from '../modules/gatk/gatk_printreads'
include {PICARD_CALCULATEHSMETRICS} from '../modules/picard/picard_calculatehsmetrics'
include {MSISENSOR2_MSI} from '../modules/msisensor2/msisensor2_msi'
include {GATK_GETSAMPLENAME} from '../modules/gatk/gatk_getsamplename'
include {GATK_MUTECT2} from '../modules/gatk/gatk_mutect2'
include {GATK_FILTERMUTECTCALLS} from '../modules/gatk/gatk_filtermutectcalls'
include {ALLELE_DEPTH_MIN_AND_AF_FROM_ADS as AD_min_AF_MUT;
         ALLELE_DEPTH_MIN_AND_AF_FROM_ADS as AD_min_AF_IND} from '../modules/utility_modules/allele_depth_min_and_AF_from_ADs'
include {SNPSIFT_ANNOTATE as ANNOTATE_AD;
         SNPSIFT_ANNOTATE as ANNOTATE_ID;
         SNPSIFT_ANNOTATE as ANNOTATE_BCF} from '../modules/snpeff_snpsift/snpsift_annotate'
include {ADD_CALLER_GATK} from '../modules/utility_modules/add_caller_gatk'
include {JOIN_ADJACENT_SNPS_AS} from '../modules/utility_modules/join_adjacent_snps_as'
include {BCF_ANNOTATE} from '../modules/bcftools/bcftools_annotate'
include {MICROINDEL_CALLING_A} from '../modules/utility_modules/microindel_calling_a'
include {MICROINDEL_CALLING_B} from '../modules/utility_modules/microindel_calling_b'
include {ADD_CALLER_PINDEL} from '../modules/utility_modules/add_caller_pindel'
include {SNPSIFT_MICROINDELS} from "${projectDir}/modules/snpeff_snpsift/snpsift_microindels"
include {SNPEFF_ANNOTATE} from "${projectDir}/modules/snpeff_snpsift/snpeff_annotate"
include {SNPSIFT_DBNSFP} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {SNPSIFT_COSMIC} from "${projectDir}/modules/snpeff_snpsift/snpsift_cosmic"
include {EXTRACT_FIELDS} from "${projectDir}/modules/snpeff_snpsift/extract_fields"
include {TMB_SCORE_PREPROCESS} from "${projectDir}/modules/utility_modules/tmb_score_preprocess"
include {TMB_SCORE} from "${projectDir}/modules/utility_modules/tmb_score"
include {SUMMARY_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_exome"
include {GATK_DEPTHOFCOVERAGE} from "${projectDir}/modules/gatk/gatk_depthofcoverage"
include {COVCALC_GATK} from "${projectDir}/modules/utility_modules/covcalc_gatk"


// log params
param_log()

// prepare reads channel
if (params.concat_lanes){
  if (params.read_type == 'PE'){
    read_ch = Channel
            .fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true, flat:true )
            .map { file, file1, file2 -> tuple(getLibraryId(file), file1, file2) }
            .groupTuple()
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}", checkExists:true, size:1 )
                .map { file, file1 -> tuple(getLibraryId(file), file1) }
                .groupTuple()
                .map{t-> [t[0], t[1].flatten()]}
  }
} else {
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
}

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}


// main workflow
workflow WES {

  // Remove `pipeline_complete.txt` from prior run, if this is a 'resume' or sample re-run. 
  // This file is used in 'on.complete' and in JAX PDX loader
  run_check = file("${params.pubdir}/pipeline_complete.txt")
  run_check.delete()

  // Create `pipeline_running.txt` used in 'on.complete' and in JAX PDX loader
  RUN_START()

  // Step 0: Concatenate Fastq files if required.
  if (params.concat_lanes){
    if (params.read_type == 'PE'){
        CONCATENATE_READS_PE(read_ch)
        read_ch = CONCATENATE_READS_PE.out.concat_fastq
    } else if (params.read_type == 'SE'){
        CONCATENATE_READS_SE(read_ch)
        read_ch = CONCATENATE_READS_SE.out.concat_fastq
    }
  }

  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)

  // Step 2: Xenome Classification A
  XENOME_CLASSIFY(QUALITY_STATISTICS.out.trimmed_fastq)

  // Step 3: Xenome Classification B
  XENOME_SORT(XENOME_CLASSIFY.out.xenome_fastq)

  // Step 4: Get Read Group Information
  READ_GROUPS(XENOME_SORT.out.sorted_fastq, "gatk")

  // Step 5: BWA-MEM Alignment
  xenome_and_rg = XENOME_SORT.out.sorted_fastq.join(READ_GROUPS.out.read_groups)
  BWA_MEM(xenome_and_rg)

  // Step 6: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.bam)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

  // Step 7: Realigner target creator and indel realigner

  dedup_and_index = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
  GATK_REALIGNERTARGETCREATOR(dedup_and_index)

  dedup_index_and_intervals = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai).join(GATK_REALIGNERTARGETCREATOR.out.intervals)
  GATK_INDELREALIGNER(dedup_index_and_intervals)

  // Step 8: Variant Pre-Processing - Part 2
  GATK_BASERECALIBRATOR(dedup_and_index)

  // Step 9: PrintReads
  dedup_index_and_grp = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai).join(GATK_BASERECALIBRATOR.out.grp)
  GATK_PRINTREADS(dedup_index_and_grp)

  // Step 10: Calculate ehsmetrics
  printreads_and_index = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai)
  PICARD_CALCULATEHSMETRICS(printreads_and_index)

  // Step 11: MSIsensor2
  MSISENSOR2_MSI(printreads_and_index)

  // Step 12: Get sample name
  GATK_GETSAMPLENAME(printreads_and_index)

  // Step 13: Mutect2
  printreads_index_and_name = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai).join(GATK_GETSAMPLENAME.out)
  GATK_MUTECT2(printreads_index_and_name)

  // Step 14 : Filter Muctect calls
  vcf_and_index = GATK_MUTECT2.out.vcf.join(GATK_MUTECT2.out.tbi)
  GATK_FILTERMUTECTCALLS(vcf_and_index)

  // Step 15 : Recompute the locus depth and Add Estimated Allele Frequency
  AD_min_AF_MUT(GATK_FILTERMUTECTCALLS.out.vcf)

  // Step 16 : Snpsift Annotate
  ANNOTATE_AD(AD_min_AF_MUT.out.vcf)

  // Step 17 : Add caller gatk
  ADD_CALLER_GATK(ANNOTATE_AD.out.vcf)

  // Step 18 : Join Adjacent SNPs AS
  join_adjacent_snps = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai).join(ADD_CALLER_GATK.out.vcf)
  JOIN_ADJACENT_SNPS_AS(join_adjacent_snps)

  // Step 19 : Bcftools Annotate
  bcf_annotate = JOIN_ADJACENT_SNPS_AS.out.vcf.join(JOIN_ADJACENT_SNPS_AS.out.tbi)
  BCF_ANNOTATE(bcf_annotate)

  // Step 20 : Snpsift Annotate
  ANNOTATE_BCF(BCF_ANNOTATE.out.vcf)

  // Step 21 : Microindel calling part 1
  microindel_call_a = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai)
  MICROINDEL_CALLING_A(microindel_call_a)

  // Step 22 : Microindel calling part 2
  MICROINDEL_CALLING_B(MICROINDEL_CALLING_A.out.vcf) 

  // Microindel filtering
  // Step 23 : Recompute the locus depth and Add Estimated Allele Frequency
  AD_min_AF_IND(MICROINDEL_CALLING_B.out.vcf)

  // Step 24 : Snpsift Annotate
  ANNOTATE_ID(AD_min_AF_IND.out.vcf)

  // Step 25 : Add caller pindel and get microIndels.DPfiltered.vcf 
  ADD_CALLER_PINDEL(ANNOTATE_ID.out.vcf)

  // Step 26 : Merge earlier annotated variants with microindels
  variants_and_microindels = ANNOTATE_BCF.out.vcf.join(ADD_CALLER_PINDEL.out.vcf)
  SNPSIFT_MICROINDELS(variants_and_microindels)

  // Variant annotation
  // Step 27: Annotation with snpsift
  SNPEFF_ANNOTATE(SNPSIFT_MICROINDELS.out.vcf)

  // Step 28: Annotation with DBNSFP
  SNPSIFT_DBNSFP(SNPEFF_ANNOTATE.out.vcf, "BOTH")

  // Step 29: Annotation with COSMIC
  SNPSIFT_COSMIC(SNPSIFT_DBNSFP.out.vcf)
  
  // Step 30: Extract required annotated fields and prepare output tables
  EXTRACT_FIELDS(SNPSIFT_COSMIC.out.vcf)

  // Step 31: Calculate TMB score from variants and microindels
  tmb_input = ADD_CALLER_GATK.out.vcf.join(ADD_CALLER_PINDEL.out.vcf)
  TMB_SCORE_PREPROCESS(tmb_input)

  tmb_input_postprocess = TMB_SCORE_PREPROCESS.out.tab.join(TMB_SCORE_PREPROCESS.out.count2).join(TMB_SCORE_PREPROCESS.out.count3)
  TMB_SCORE(tmb_input_postprocess)

  // Step 32: Summary statistics

  fq_alignment_metrics = QUALITY_STATISTICS.out.quality_stats.join(PICARD_MARKDUPLICATES.out.dedup_metrics).join(PICARD_CALCULATEHSMETRICS.out.hsmetrics)
  SUMMARY_STATS(fq_alignment_metrics)

  depth_of_coverage_hex = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai)
  GATK_DEPTHOFCOVERAGE(depth_of_coverage_hex, params.hex_genes)
  COVCALC_GATK(GATK_DEPTHOFCOVERAGE.out.txt, "HEX")
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