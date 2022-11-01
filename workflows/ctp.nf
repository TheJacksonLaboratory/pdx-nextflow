#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {XENOME_CLASSIFY} from   "${projectDir}/modules/xenome/xenome"
include {FASTQ_SORT as XENOME_SORT} from   "${projectDir}/modules/fastq-tools/fastq-tools_sort"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {GATK_REALIGNERTARGETCREATOR} from "${projectDir}/modules/gatk/gatk_realignertargetcreator"
include {GATK_INDELREALIGNER} from "${projectDir}/modules/gatk/gatk_indelrealigner"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_PRINTREADS} from "${projectDir}/modules/gatk/gatk_printreads"
include {PICARD_CALCULATEHSMETRICS} from "${projectDir}/modules/picard/picard_calculatehsmetrics"

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
workflow CTP {

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

  // Step 2: Xenome Classification
  XENOME_CLASSIFY(QUALITY_STATISTICS.out.trimmed_fastq)

  // Step 3: Sorting Post-Xenome Reads
  XENOME_SORT(XENOME_CLASSIFY.out.xenome_fastq)

  // Step 4: Get read group information
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

  // Step 10: Calculate depth metrics
  printreads_and_index = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai)
  PICARD_CALCULATEHSMETRICS(printreads_and_index)
}