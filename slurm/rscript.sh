#!/usr/bin/env bash

if [[ -z "$DOCKER_IMG" ]]
then
    echo "Error: DOCKER_IMG environment variable empty" 1>&2
    exit 1
fi

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
