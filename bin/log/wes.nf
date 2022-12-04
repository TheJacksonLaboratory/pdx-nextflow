def param_log(){
log.info """
______________________________________________________

        WHOLE EXOME SEQUENCING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--ref_prefix                    ${params.ref_prefix}
--gold_std_indels               ${params.gold_std_indels}
--known_indels                  ${params.known_indels}
--targets_gatk                  ${params.targets_gatk}
--dbsnp                         ${params.dbsnp}
--targets_picard                ${params.targets_picard}
--targets_exons_picard          ${params.targets_exons_picard}
--hex_percentage                ${params.hex_percentage}
--msisensor_model               ${params.msisensor_model}
--exac_ref                      ${params.exac_ref}
--samp_ploidy                   ${params.samp_ploidy}
--minDP                         ${params.minDP}
--fa2bit                        ${params.fa2bit}
--snpEff_config                 ${params.snpEff_config}
--snpEff_data                   ${params.snpEff_data}
--hgvs_data                     ${params.hgvs_data}
--ensembl_transcript            ${params.ensembl_transcript}
--dbNSFP                        ${params.dbNSFP}
--Cosmic_newer                  ${params.Cosmic_newer}
--Cosmic_older                  ${params.Cosmic_older}
--rec_var                       ${params.rec_var}
--bins_hexcoverage              ${params.bins_hexcoverage}
--hex_genes                     ${params.hex_genes}
--preserve_work                 ${params.preserve_work}

Project Directory: ${projectDir}
______________________________________________________
"""
}