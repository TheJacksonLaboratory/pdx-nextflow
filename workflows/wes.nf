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
         SNPSIFT_ANNOTATE as ANNOTATE_BCF} from '../modules/snpeff_snpsift/snpsift_annotate'
include {ADD_CALLER_GATK} from '../modules/utility_modules/add_caller_gatk'
include {JOIN_ADJACENT_SNPS_AS} from '../modules/utility_modules/join_adjacent_snps_as'
include {BCF_ANNOTATE} from '../modules/bcftools/bcftools_annotate'


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

  // Step 6: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.bam)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

  // Step 7: Realigner target creator and indel realigner
  GATK_REALIGNERTARGETCREATOR(PICARD_MARKDUPLICATES.out.dedup_bam, PICARD_MARKDUPLICATES.out.dedup_bai)
  GATK_INDELREALIGNER(PICARD_MARKDUPLICATES.out.dedup_bam, PICARD_MARKDUPLICATES.out.dedup_bai,
                      GATK_REALIGNERTARGETCREATOR.out.intervals)

  // Step 8: Variant Pre-Processing - Part 2
  GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam, PICARD_MARKDUPLICATES.out.dedup_bai)

  // Step 9: PrintReads
  GATK_PRINTREADS(PICARD_MARKDUPLICATES.out.dedup_bam, PICARD_MARKDUPLICATES.out.dedup_bai, GATK_BASERECALIBRATOR.out.grp)

  // Step 10: Calculate ehsmetrics
  PICARD_CALCULATEHSMETRICS(GATK_PRINTREADS.out.bam, GATK_PRINTREADS.out.bai)

  // Step 11: MSIsensor2
  MSISENSOR2_MSI(GATK_PRINTREADS.out.bam, GATK_PRINTREADS.out.bai)

  // Step 12: Get sample name
  GATK_GETSAMPLENAME(GATK_PRINTREADS.out.bam, GATK_PRINTREADS.out.bai)

  // Step 13: Mutect2
  GATK_MUTECT2(GATK_PRINTREADS.out.bam, GATK_PRINTREADS.out.bai, GATK_GETSAMPLENAME.out)

  // Step 14 : Filter Muctect calls
  GATK_FILTERMUTECTCALLS(GATK_MUTECT2.out.vcf, GATK_MUTECT2.out.tbi)

  // Step 15 : Recompute the locus depth and Add Estimated Allele Frequency
  AD_min_AF_MUT(GATK_FILTERMUTECTCALLS.out.vcf)

  // Step 16 : Snpsift Annotate
  ANNOTATE_AD(AD_min_AF_MUT.out.vcf)

  // Step 17 : Add caller gatk
  ADD_CALLER_GATK(ANNOTATE_AD.out.vcf)

  // Step 18 : Join Adjacent SNPs AS
  JOIN_ADJACENT_SNPS_AS(GATK_PRINTREADS.out.bam, GATK_PRINTREADS.out.bai, ADD_CALLER_GATK.out.vcf)

  // Step 19 : Bcftools Annotate
  BCF_ANNOTATE(JOIN_ADJACENT_SNPS_AS.out.vcf, JOIN_ADJACENT_SNPS_AS.out.tbi)

  // Step 20 : Snpsift Annotate
  ANNOTATE_BCF(BCF_ANNOTATE.out.vcf)



}

