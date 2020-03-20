#!/usr/bin/env bash
source config.sh
export ENTREZ_KEY

# run data generation steps with each waiting for the previous to finish
STEP1_JOB_ID=$(sbatch --job-name="ddh-step1" run-slurm-step.sh "1")
STEP2_JOB_ID=$(sbatch --job-name="ddh-step2" --dependency=afterany:$STEP1_JOB_ID run-slurm-step.sh --mem=180G "2")
STEP3_JOB_ID=$(sbatch --job-name="ddh-step3" --dependency=afterany:$STEP2_JOB_ID --cpus-per-task=10 run-slurm-step.sh "3")
STEP4_JOB_ID=$(sbatch --job-name="ddh-step4" --dependency=afterany:$STEP3_JOB_ID run-slurm-step.sh "4")

# upload results when last data generation step finishes
sbatch --job-name="ddh-upload" --dependency=afterany:$STEP4_JOB_ID upload-slurm.sh
