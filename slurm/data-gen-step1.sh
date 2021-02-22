#!/usr/bin/env bash
#SBATCH --job-name=ddh-1
#SBATCH --mem=64G
#SBATCH --output=logs/ddh-1-%j.out

source config.sh

DDH_STEP=1 slurm/rscript.sh code/generate_data_private.R
