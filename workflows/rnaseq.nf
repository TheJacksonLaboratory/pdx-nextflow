#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from '../bin/shared/getLibraryId.nf'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {QUALITY_STATISTICS} from '../modules/utility_modules/quality_stats'
include {CLASSIFICATION_A} from   '../modules/xenome/xenome_classification_a'
include {CLASSIFICATION_B} from   '../modules/fastq-tools/fastq-tools_sort'
include {RSEM_ALIGNMENT_EXPRESSION} from '../modules/rsem/rsem_alignment_expression'
include {ADD_GENE_NAME_NORM} from '../modules/perl/perl_add_gene_name_and_normalization'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {PICARD_ADDORREPLACEREADGROUPS} from '../modules/picard/picard_addorreplacereadgroups'
include {PICARD_REORDERSAM} from '../modules/picard/picard_reordersam'
include {PICARD_COLLECTRNASEQMETRICS} from '../modules/picard/picard_collectrnaseqmetrics'
include {PICARD_SORTSAM} from '../modules/picard/picard_sortsam'
include {GATK_DEPTHOFCOVERAGE as GATK_DEPTHOFCOVERAGE_CTP;
         GATK_DEPTHOFCOVERAGE as GATK_DEPTHOFCOVERAGE_PROBES} from '../modules/gatk/gatk_depthofcoverage'
include {FORMAT_GATK as FORMAT_GATK_CTP;
         FORMAT_GATK as FORMAT_GATK_PROBES} from '../modules/utility_modules/rna_format_gatk'
include {COVCALC_GATK as COVCALC_GATK_CTP;
         COVCALC_GATK as COVCALC_GATK_PROBES} from '../modules/utility_modules/rna_covcalc_gatk'
include {CLASSIFIER_COVERAGE} from '../modules/utility_modules/rna_classifier_coverage'
include {RNA_SUMMARY_STATS} from '../modules/utility_modules/aggregate_stats_rna'


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
workflow RNASEQ {

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
  CLASSIFICATION_A(QUALITY_STATISTICS.out.trimmed_fastq)

  // Step 3: Xenome Classification B
  CLASSIFICATION_B(CLASSIFICATION_A.out.xenome_fastq)

  // Step 4: RSEM
  RSEM_ALIGNMENT_EXPRESSION(CLASSIFICATION_B.out.sorted_fastq)

  // Step 5: JOIN
  data = RSEM_ALIGNMENT_EXPRESSION.out.rsem_genes.join(RSEM_ALIGNMENT_EXPRESSION.out.rsem_isoforms)

  // Step 6: ADD GENE NAME AND NORMALIZATION
  ADD_GENE_NAME_NORM(data)

  // Step 7: Get Read Group Information
  READ_GROUPS(read_ch, "picard")


  // Step 8: Picard Alignment Metrics
  PICARD_ADDORREPLACEREADGROUPS(READ_GROUPS.out.read_groups,
                                RSEM_ALIGNMENT_EXPRESSION.out.bam)

  // Step 9: Picard Alignment Metrics
  PICARD_REORDERSAM(PICARD_ADDORREPLACEREADGROUPS.out.bam)

  // Step 10: Picard Alignment Metrics
  PICARD_SORTSAM(PICARD_REORDERSAM.out.bam)

  // Step 11: Picard Alignment Metrics
  PICARD_COLLECTRNASEQMETRICS(PICARD_SORTSAM.out.bam)  

  // Step 12: GATK Coverage Stats (CTP)
  GATK_DEPTHOFCOVERAGE_CTP(PICARD_SORTSAM.out.bam, PICARD_SORTSAM.out.bai, params.ctp_genes)
  FORMAT_GATK_CTP(GATK_DEPTHOFCOVERAGE_CTP.out.txt, params.ctp_genes)
  COVCALC_GATK_CTP(FORMAT_GATK_CTP.out.txt, "CTP")

  // Step 13: GATK Coverage Stats (PROBES)
  GATK_DEPTHOFCOVERAGE_PROBES(PICARD_SORTSAM.out.bam, PICARD_SORTSAM.out.bai, params.probes)
  FORMAT_GATK_PROBES(GATK_DEPTHOFCOVERAGE_PROBES.out.txt, params.probes)
  COVCALC_GATK_PROBES(FORMAT_GATK_PROBES.out.txt, "PROBES")

  // Step 14: Classifier and Coverage
  CLASSIFIER_COVERAGE(ADD_GENE_NAME_NORM.out.norm_gene_results)

  // Step 15: Summary Stats
  RNA_SUMMARY_STATS(RSEM_ALIGNMENT_EXPRESSION.out.rsem_stats,
                    QUALITY_STATISTICS.out.quality_stats,
                    CLASSIFICATION_A.out.xenome_stats,
                    PICARD_COLLECTRNASEQMETRICS.out.picard_metrics)


}
