#!/usr/bin/env bash

# setup DDH_EMAIL environment variable
source config.sh

# verify we have config.sh setup correctly
if [[ -z "$DOCKER_IMG" ]]
then
    echo "Error: DOCKER_IMG environment variable empty" 1>&2
    exit 1
fi
if [[ -z "$DDS_PROJECT" ]]
then
    echo "Error: DDS_PROJECT environment variable empty" 1>&2
    exit 1
fi

# create flags that will send an email when a job fails or when a job ends(fail or succeed)
FAIL_EMAIL_FLAGS="--mail-type=FAIL --mail-user=$DDH_EMAIL"
DONE_EMAIL_FLAGS="--mail-type=END --mail-user=$DDH_EMAIL"

# make sure we have data and logs directories created
mkdir -p data
mkdir -p logs

# submit sbatch jobs that wait on the previous  job to finish
JOB_ID=""
for SCRIPT in data-gen-step1.sh \
              data-gen-step2.sh \
              data-gen-step3.sh \
              data-gen-step4-pos-subsets.sh \
              data-gen-step4-neg-subsets.sh \
              data-gen-step4-merge-subsets.sh
do
    JOB_FLAGS=${FAIL_EMAIL_FLAGS}
    if [ "$JOB_ID" != "" ]
    then
      JOB_FLAGS="$FAIL_EMAIL_FLAGS --dependency=afterok:$JOB_ID"
    fi
    echo "sbatch ${JOB_FLAGS} --parsable slurm/${SCRIPT}"
    JOB_ID=$(sbatch ${JOB_FLAGS} --parsable slurm/${SCRIPT})
    echo "Created job $JOB_ID for $SCRIPT"
done

echo "sbatch ${DONE_EMAIL_FLAGS} --job-name="ddh-upload" --dependency=afterok:$JOB_ID slurm/upload-data.sh"
UPLOAD_JOB_ID=$(sbatch ${DONE_EMAIL_FLAGS} --job-name="ddh-upload" --dependency=afterok:$JOB_ID slurm/upload-data.sh)
echo "Created job $JOB_ID for slurm/upload-data.sh"
