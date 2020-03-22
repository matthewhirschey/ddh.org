#!/usr/bin/env bash
source config.sh
export ENTREZ_KEY

STEP1_FLAGS="--mem=64G"
STEP2_FLAGS="--mem=212G"
STEP3_FLAGS="--mem=32G"
STEP4_FLAGS="--mem=32G --cpus-per-task=10"
STD_FLAGS="--parsable"
if [ -v "$DDH_EMAIL" ]
then
  STD_FLAGS="$STD_FLAGS --mail-type=DONE --mail-user=$DDH_EMAIL"
fi

# run data generation steps with each waiting for the previous to finish
STEP1_JOB_ID=$(sbatch $STD_FLAGS --job-name="ddh-gens1" $STEP1_FLAGS run-slurm-step.sh "gens")
STEP2_JOB_ID=$(sbatch $STD_FLAGS --job-name="ddh-pubm2" --dependency=afterok:$STEP1_JOB_ID $STEP2_FLAGS run-slurm-step.sh "pubm")
STEP3_JOB_ID=$(sbatch $STD_FLAGS --job-name="ddh-tbls3" --dependency=afterok:$STEP2_JOB_ID $STEP3_FLAGS run-slurm-step.sh "tbls")
STEP4_JOB_ID=$(sbatch $STD_FLAGS --job-name="ddh-path4" --dependency=afterok:$STEP3_JOB_ID $STEP4_FLAGS run-slurm-step.sh "path")

if [ -v "$DDH_UPLOAD_RESULTS" ]
then
   # upload results when last data generation step finishes
   sbatch $STD_FLAGS --job-name="ddh-upload" --dependency=afterok:$STEP4_JOB_ID upload-slurm.sh
fi
