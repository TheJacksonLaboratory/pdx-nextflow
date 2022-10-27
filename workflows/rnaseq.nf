#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {XENOME_CLASSIFY} from   "${projectDir}/modules/xenome/xenome"
include {FASTQ_SORT as XENOME_SORT} from   "${projectDir}/modules/fastq-tools/fastq-tools_sort"
include {RSEM_ALIGNMENT_EXPRESSION} from "${projectDir}/modules/rsem/rsem_alignment_expression"
include {ADD_GENE_NAME_NORM} from "${projectDir}/modules/perl/perl_add_gene_name_and_normalization"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {PICARD_ADDORREPLACEREADGROUPS} from "${projectDir}/modules/picard/picard_addorreplacereadgroups"
include {PICARD_REORDERSAM} from "${projectDir}/modules/picard/picard_reordersam"
include {PICARD_COLLECTRNASEQMETRICS} from "${projectDir}/modules/picard/picard_collectrnaseqmetrics"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {GATK_DEPTHOFCOVERAGE as GATK_DEPTHOFCOVERAGE_CTP;
         GATK_DEPTHOFCOVERAGE as GATK_DEPTHOFCOVERAGE_PROBES} from "${projectDir}/modules/gatk/gatk_depthofcoverage"
include {FORMAT_GATK as FORMAT_GATK_CTP;
         FORMAT_GATK as FORMAT_GATK_PROBES} from "${projectDir}/modules/utility_modules/rna_format_gatk"
include {COVCALC_GATK as COVCALC_GATK_CTP;
         COVCALC_GATK as COVCALC_GATK_PROBES} from "${projectDir}/modules/utility_modules/rna_covcalc_gatk"
include {CLASSIFIER_COVERAGE} from "${projectDir}/modules/utility_modules/rna_classifier_coverage"
include {RNA_SUMMARY_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_rna"


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
  XENOME_CLASSIFY(QUALITY_STATISTICS.out.trimmed_fastq)

  // Step 3: Xenome Classification B
  XENOME_SORT(XENOME_CLASSIFY.out.xenome_fastq)

  // Step 4: RSEM
  RSEM_ALIGNMENT_EXPRESSION(XENOME_SORT.out.sorted_fastq)

  // Step 5: JOIN
  data = RSEM_ALIGNMENT_EXPRESSION.out.rsem_genes.join(RSEM_ALIGNMENT_EXPRESSION.out.rsem_isoforms)

  // Step 6: ADD GENE NAME AND NORMALIZATION
  ADD_GENE_NAME_NORM(data)

  // Step 7: Get Read Group Information
  READ_GROUPS(read_ch, "picard")


  // Step 8: Picard Alignment Metrics
  add_replace_groups = READ_GROUPS.out.read_groups.join(RSEM_ALIGNMENT_EXPRESSION.out.bam)
  PICARD_ADDORREPLACEREADGROUPS(add_replace_groups)

  // Step 9: Picard Alignment Metrics
  PICARD_REORDERSAM(PICARD_ADDORREPLACEREADGROUPS.out.bam)

  // Step 10: Picard Alignment Metrics
  PICARD_SORTSAM(PICARD_REORDERSAM.out.bam)

  // Step 11: Picard Alignment Metrics
  PICARD_COLLECTRNASEQMETRICS(PICARD_SORTSAM.out.bam)  

  // Step 12: GATK Coverage Stats (CTP)
  depth_of_coverage_ctp = PICARD_SORTSAM.out.bam.join(PICARD_SORTSAM.out.bai)
  GATK_DEPTHOFCOVERAGE_CTP(depth_of_coverage_ctp, params.ctp_genes)
  FORMAT_GATK_CTP(GATK_DEPTHOFCOVERAGE_CTP.out.txt, params.ctp_genes)
  COVCALC_GATK_CTP(FORMAT_GATK_CTP.out.txt, "CTP")

  // Step 13: GATK Coverage Stats (PROBES)
  depth_of_coverage_probes = PICARD_SORTSAM.out.bam.join(PICARD_SORTSAM.out.bai)
  GATK_DEPTHOFCOVERAGE_PROBES(depth_of_coverage_probes, params.probes)
  FORMAT_GATK_PROBES(GATK_DEPTHOFCOVERAGE_PROBES.out.txt, params.probes)
  COVCALC_GATK_PROBES(FORMAT_GATK_PROBES.out.txt, "PROBES")

  // Step 14: Classifier and Coverage
  CLASSIFIER_COVERAGE(ADD_GENE_NAME_NORM.out.norm_gene_results)

  // Step 15: Summary Stats
  aggregate_stats_rna = RSEM_ALIGNMENT_EXPRESSION.out.rsem_stats.join(QUALITY_STATISTICS.out.quality_stats).join(XENOME_CLASSIFY.out.xenome_stats).join(PICARD_COLLECTRNASEQMETRICS.out.picard_metrics)
  RNA_SUMMARY_STATS(aggregate_stats_rna)


}
