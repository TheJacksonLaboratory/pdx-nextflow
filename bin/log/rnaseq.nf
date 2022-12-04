def param_log(){
log.info """
______________________________________________________

                RNA-seq PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--read_type                     ${params.read_type}
--read_prep                     ${params.read_prep}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--seed_length                   ${params.seed_length}
--rsem_aligner                  ${params.rsem_aligner}
--hsa_accession                 ${params.hsa_accession}
--ref_prefix                    ${params.ref_prefix}
--ref_fa                        ${params.ref_fa}
--ref_flat                      ${params.ref_flat}
--ribo_intervals                ${params.ribo_intervals}
--probes                        ${params.probes}
--ctp_genes                     ${params.ctp_genes}
--classifier_table              ${params.classifier_table}
--preserve_work                 ${params.preserve_work}

Project Directory: ${projectDir}
______________________________________________________
"""
}