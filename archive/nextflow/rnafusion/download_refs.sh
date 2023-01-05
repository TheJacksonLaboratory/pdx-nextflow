#!/bin/bash

#SBATCH --job-name=SlurmJob
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=mike.lloyd@jax.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH --output=slurm-%j.out

nextflow run ./download-references.nf -profile singularity --download_all --cosmic_usr mike.lloyd@jax.org --cosmic_passwd YSYLTvNy72fvxg!
