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
include {MSISENSOR2_MSI} from "${projectDir}/modules/msisensor2/msisensor2_msi"
include {GATK_GETSAMPLENAME} from "${projectDir}/modules/gatk/gatk_getsamplename"
include {GATK_MUTECT2_CTP} from "${projectDir}/modules/gatk/gatk_mutect2_ctp"
include {GATK_FILTERMUTECTCALLS_CTP} from "${projectDir}/modules/gatk/gatk_filtermutectcalls_ctp"
include {ALLELE_DEPTH_MIN_AND_AF_FROM_ADS as AD_MIN_AF_MUT;
         ALLELE_DEPTH_MIN_AND_AF_FROM_ADS as AD_MIN_AF_IND} from "${projectDir}/modules/utility_modules/allele_depth_min_and_AF_from_ADs"
include {SNPSIFT_ANNOTATE as ANNOTATE_AD;
         SNPSIFT_ANNOTATE as ANNOTATE_ID;
         SNPSIFT_ANNOTATE as ANNOTATE_BCF} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {ADD_CALLER_GATK} from "${projectDir}/modules/utility_modules/add_caller_gatk"
include {JOIN_ADJACENT_SNPS_AS} from "${projectDir}/modules/utility_modules/join_adjacent_snps_as"
include {BCF_ANNOTATE} from "${projectDir}/modules/bcftools/bcftools_annotate"


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

  // Step 11: MSIsensor2
  MSISENSOR2_MSI(printreads_and_index)

  // Step 12: Get sample name
  GATK_GETSAMPLENAME(printreads_and_index)

  // Step 13: Mutect2
  printreads_index_and_name = GATK_PRINTREADS.out.bam.join(GATK_PRINTREADS.out.bai).join(GATK_GETSAMPLENAME.out)
  GATK_MUTECT2_CTP(printreads_index_and_name)

  // Step 14 : Filter Mutect calls
  vcf_and_index = GATK_MUTECT2_CTP.out.vcf.join(GATK_MUTECT2_CTP.out.tbi)
  GATK_FILTERMUTECTCALLS_CTP(vcf_and_index)

  // Step 15 : Recompute the locus depth and Add Estimated Allele Frequency
  AD_MIN_AF_MUT(GATK_FILTERMUTECTCALLS_CTP.out.vcf)

  // Step 16 : Snpsift Annotate
  ANNOTATE_AD(AD_MIN_AF_MUT.out.vcf)

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

}