#!/usr/bin/env bash
DOCKER_IMG="docker://matthewhirschey/ddh:latest"

# numeric step to run
export DDH_STEP=$1

# make sure we have a data directory created
mkdir data

# setup singularity environment variables 
SINGULARITY_BASE="singularity"
export SINGULARITY_TMPDIR="${SINGULARITY_BASE}/tmp"
export SINGULARITY_CACHEDIR="${SINGULARITY_BASE}/cache"
export SINGULARITY_IMAGEDIR="${SINGULARITY_BASE}/images"
mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $SINGULARITY_IMAGEDIR

# pull singularity image if it doesn't exist
if [ ! -f "$SINGULARITY_IMAGEDIR/ddh.sif" ]
then
  singularity pull $SINGULARITY_IMAGEDIR/ddh.sif $DOCKER_IMG
fi

# run data generation step
singularity exec $SINGULARITY_IMAGEDIR/ddh.sif Rscript code/generate_data.R
