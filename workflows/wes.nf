#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from '../bin/shared/getLibraryId.nf'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {QUALITY_STATISTICS} from '../modules/utility_modules/quality_stats'
include {XENOME_CLASSIFY} from   '../modules/xenome/xenome'
include {FASTQ_SORT as XENOME_SORT} from   '../modules/fastq-tools/fastq-tools_sort'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {BWA_MEM} from '../modules/bwa/bwa_mem'

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
  BWA_MEM(XENOME_SORT.out.sorted_fastq, READ_GROUPS.out.read_groups )

}

