#!/usr/bin/env bash

DOCKER_IMG="docker://matthewhirschey/ddh.com:latest"

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

singularity exec $SINGULARITY_IMAGEDIR/ddh.sif Rscript "$@"
