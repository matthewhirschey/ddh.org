#!/usr/bin/env bash
#SBATCH --mem=32G

source config.sh

make container_image
RSCRIPT_CMD="singularity exec singularity/depmap.sif Rscript" make
