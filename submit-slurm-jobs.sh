#!/usr/bin/env bash
source config.sh
export ENTREZ_KEY

STEP1_FLAGS="--mem=64G"
STEP2_FLAGS="--mem=212G"
STEP3_FLAGS="--mem=32G"
STEP4_FLAGS="--mem=32G --cpus-per-task=10"

# run data generation steps with each waiting for the previous to finish
STEP1_JOB_ID=$(sbatch --parsable --job-name="ddh-step1" $STEP1_FLAGS run-slurm-step.sh "1")
STEP2_JOB_ID=$(sbatch --parsable --job-name="ddh-step2" --dependency=afterok:$STEP1_JOB_ID $STEP2_FLAGS run-slurm-step.sh "2")
STEP3_JOB_ID=$(sbatch --parsable --job-name="ddh-step3" --dependency=afterok:$STEP2_JOB_ID $STEP3_FLAGS run-slurm-step.sh "3")
STEP4_JOB_ID=$(sbatch --parsable --job-name="ddh-step4" --dependency=afterok:$STEP3_JOB_ID $STEP4_FLAGS run-slurm-step.sh "4")

# upload results when last data generation step finishes
sbatch --job-name="ddh-upload" --dependency=afterok:$STEP4_JOB_ID upload-slurm.sh
