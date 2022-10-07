#!/usr/bin/env bash
#SBATCH --job-name=ddh-5
#SBATCH --mem=112G
#SBATCH --output=logs/ddh-5-%j.out

source config.sh
DDH_STEP=5 slurm/rscript.sh code/generate_methods.R
