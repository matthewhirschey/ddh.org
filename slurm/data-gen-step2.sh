#!/usr/bin/env bash
#SBATCH --job-name=ddh-2
#SBATCH --mem=270G
#SBATCH --output=logs/ddh-2-%j.out

source config.sh
DDH_STEP=2 slurm/rscript.sh code/generate_data.R
