#!/usr/bin/env bash
#SBATCH --job-name=ddh-4-merge
#SBATCH --mem=32G
#SBATCH --output=logs/ddh-4-merge-%j.out

source config.sh

DDH_STEP="4.merge" slurm/rscript.sh code/generate_data_private.R
