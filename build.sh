#!/usr/bin/env bash
#SBATCH --mem=32G

source config.sh

RSCRIPT_CMD="singularity exec singularity/depmap.sif Rscript" make gene_summary
