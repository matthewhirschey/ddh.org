#!/usr/bin/env bash
#SBATCH --job-name=ddh-3
#SBATCH --mem=96G
#SBATCH --output=logs/ddh-3-%j.out

source config.sh

export ENTREZ_KEY
DDH_STEP=3 slurm/rscript.sh code/generate_data.R
