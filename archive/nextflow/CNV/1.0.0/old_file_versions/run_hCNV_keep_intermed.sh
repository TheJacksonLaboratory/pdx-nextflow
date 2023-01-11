#!/bin/bash
#SBATCH --job-name=nf_hCNV
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@jax.org
#SBATCH -q batch
#SBATCH -t 24:00:00
#SBATCH --mem=2000



~/bin/nextflow \
-c ./params.config \
run \
/pdx/pdx_resource_service/elion/pipelines/CNV/1.0.0/hCNV_keep_intermed.nf \
-profile slurm,singularity
