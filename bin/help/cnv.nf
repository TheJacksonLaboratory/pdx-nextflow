def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow -c path/to/params.cfg run path/to/pipeline.nf -profile slurm, singularity
            (The params.cfg file needs to have the following mandatory parameters
             OR they need to specified on the command line.)

    Mandatory:
        --celInput                absolute path to input CEL.
        --hapmap_dat              directory containing hapmap data from NCBI
        --hapmap_fm               file containing hapmap data for females
        --hapmap_m                file containing hapmap data for males
        --genoclust               genotype clustering file
        --snp6chip                file containing description of the probe sets on the chip
        --snp6chip_birdseed_mod   file containing birdseed models
        --snp6chip_specsnps       file containing chromosme X (non-pseudo-autosomal), chromosome Y, and mitochondrial SNPs
        --hapmap_norm_target      file containing hapmap sample data used for normalization
        --gw6_pfb_file            file containing annotated marker positions
        --SNPpos                  file containing annotated marker positions
        --GC                      file containing GC% for annotated marker positions
        --chr_arm                 file containing positions of chromosome arms
        --exp_mart_genes          file containing genes exported to PDX data mart database

    Optional:
        --preserve_work         if 'yes' then don't delete the ./work dir at end

    """.stripIndent()
}
