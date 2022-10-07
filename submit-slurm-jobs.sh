#!/usr/bin/env bash
# sbatch script that runs the data generation pipeline
# There is an optional positional argument for the restart step.
# Usage:
# sbatch submit-slurm-jobs.sh [restartstep]
#     where restartstep is 1, 2, 3, 4, or 5, default is 1 which runs all steps
#SBATCH --output=logs/submit-slurm-jobs-%j.out

# fetch restart step from positional argument - default to 1
RESTART_STEP=${1:-1}

# stop if a command fails (non-zero exit status)
set -e

echo "Starting DDH data generation"

# setup DDH_EMAIL environment variable
echo "Setup: Loading config.sh"
echo "    source config.sh"
source config.sh

# Default PULL_IMAGE to "Y"
if [[ -z "$PULL_IMAGE" ]]
then
   export PULL_IMAGE=Y
fi

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

if [ "$PULL_IMAGE" == "Y" ]
then
    # setup singularity environment variables
    SINGULARITY_BASE="$(pwd)/singularity"
    export SINGULARITY_TMPDIR="${SINGULARITY_BASE}/tmp"
    export SINGULARITY_CACHEDIR="${SINGULARITY_BASE}/cache"
    export SINGULARITY_IMAGEDIR="${SINGULARITY_BASE}/images"
    mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $SINGULARITY_IMAGEDIR

    # delete the singularity image if it exists
    if [ -f "$SINGULARITY_IMAGEDIR/ddh.sif" ]
    then
      rm $SINGULARITY_IMAGEDIR/ddh.sif
    fi

    # pull the image
    singularity pull $SINGULARITY_IMAGEDIR/ddh.sif $DOCKER_IMG
fi

# make sure we have data and logs directories created
echo "Setup: Creating data and logs directories if they do not exist."
mkdir -p data
mkdir -p logs

if [ "$RESTART_STEP" -le 1 ]
then
    echo "Running step 1"
    echo "    sbatch slurm/data-gen-step1.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step1.sh
fi

if [ "$RESTART_STEP" -le 2 ]
then
    echo "Running step 2"
    echo "    sbatch slurm/data-gen-step2.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step2.sh
fi

if [ "$RESTART_STEP" -le 3 ]
then
    echo "Running step 3"
    echo "    sbatch slurm/data-gen-step3.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step3.sh
fi

if [ "$RESTART_STEP" -le 4 ]
then
    echo "Running step 4"
    echo "    sbatch slurm/data-gen-step4-pos-subsets.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step4-pos-subsets.sh

    echo "    sbatch slurm/data-gen-step4-neg-subsets.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step4-neg-subsets.sh

    echo "    sbatch slurm/data-gen-step4-merge-subsets.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step4-merge-subsets.sh
fi

if [ "$RESTART_STEP" -le 5 ]
then
    echo "Running step 5 - generate methods"
    echo "    sbatch slurm/data-gen-step5.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step5.sh
fi

if [ "$RESTART_STEP" -le 6 ]
then
    echo "Running step 6"
    echo "    sbatch slurm/data-gen-step6.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step6.sh
fi

if [ "$RESTART_STEP" -le 7 ]
then
    echo "Running step 7 - generate cards"
    echo "    sbatch slurm/data-gen-step7.sh"
    sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-step7.sh
fi

echo "    sbatch slurm/data-gen-create-test-data.sh"
sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/data-gen-create-test-data.sh

echo "Zipping images"
echo "    sbatch slurm/zip-images.sh"
sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/zip-images.sh

echo "Uploading results"
echo "    sbatch slurm/upload-data.sh"
sbatch --wait --mail-type=END --mail-user=$DDH_EMAIL slurm/upload-data.sh

echo "Finished DDH data generation"
