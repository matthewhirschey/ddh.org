#!/usr/bin/env bash
#SBATCH --job-name=ddh-4-neg
#SBATCH --mem=32G
#SBATCH --array=1-200%5
#SBATCH --output=logs/ddh-4-neg-%A_%a.out

source config.sh
DDH_STEP="4.neg" DDH_IDX=$SLURM_ARRAY_TASK_ID slurm/rscript.sh code/generate_data.R
