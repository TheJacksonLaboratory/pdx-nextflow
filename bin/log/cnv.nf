
def param_log(){
log.info """
______________________________________________________

                CNV PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--celInput                      ${params.celInput}
--hapmap_dat                    ${params.hapmap_dat}
--hapmap_fm                     ${params.hapmap_fm}
--hapmap_m                      ${params.hapmap_m}
--genoclust                     ${params.genoclust}
--snp6chip                      ${params.snp6chip}
--snp6chip_birdseed_mod         ${params.snp6chip_birdseed_mod}
--snp6chip_specsnps             ${params.snp6chip_specsnps}
--hapmap_norm_target            ${params.hapmap_norm_target}
--gw6_pfb_file                  ${params.gw6_pfb_file}
--SNPpos                        ${params.SNPpos}
--GC                            ${params.GC}
--chr_arm                       ${params.chr_arm}
--exp_mart_genes                ${params.exp_mart_genes}
--preserve_work                 ${params.preserve_work}


Project Directory: ${projectDir}
______________________________________________________
"""
}