#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {param_log} from "${projectDir}/bin/log/rnafusion.nf"
include {RUN_START} from "${projectDir}/bin/shared/run_start"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {GUNZIP} from "${projectDir}/modules/utility_modules/gunzip"
include {XENOME_CLASSIFY} from   "${projectDir}/modules/xenome/xenome"
include {FASTQ_SORT as XENOME_SORT} from   "${projectDir}/modules/fastq-tools/fastq-tools_sort"
include {STAR_FUSION as STAR_FUSION} from "${projectDir}/modules/star-fusion/star-fusion"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {FUSION_REPORT} from "${projectDir}/modules/fusion-report/fusion_report"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

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
workflow RNAFUSION {

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

  // Step 00: Unzip Reads if needed.
  GUNZIP(read_ch)

  // Step 1: Xenome Classification A
  XENOME_CLASSIFY(GUNZIP.out.gunzip_fastq)

  // Step 2: Xenome Classification B
  XENOME_SORT(XENOME_CLASSIFY.out.xenome_fastq)

  // Step 3: Star-fusion
  STAR_FUSION(XENOME_SORT.out.sorted_fastq)

  // Step 4: FASTQC
  FASTQC(read_ch)

  // Step 5: Fusion Reporter
  FUSION_REPORT(STAR_FUSION.out.star_fusion_fusions)

  // Step 6: MultiQC
  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(FUSION_REPORT.out.summary_fusions_mq.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )
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