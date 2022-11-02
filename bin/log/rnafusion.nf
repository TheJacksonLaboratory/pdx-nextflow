def param_log(){
log.info """
______________________________________________________

                RNA FUSION PARAMETER LOG

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
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_prefix                    ${params.ref_prefix}
--star_fusion_ref               ${params.star_fusion_ref}
--star_index                    ${params.star_index}
--star_fusion_opt               ${params.star_fusion_opt}

Project Directory: ${projectDir}
______________________________________________________
"""
}