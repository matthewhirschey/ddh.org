#!/usr/bin/env bash

# stop if a command fails (non-zero exit status)
set -e

if [[ -z "$DOCKER_IMG" ]]
then
    echo "Error: DOCKER_IMG environment variable empty" 1>&2
    exit 1
fi

# setup singularity environment variables
SINGULARITY_BASE="$(pwd)/singularity"
export SINGULARITY_TMPDIR="${SINGULARITY_BASE}/tmp"
export SINGULARITY_CACHEDIR="${SINGULARITY_BASE}/cache"
export SINGULARITY_IMAGEDIR="${SINGULARITY_BASE}/images"
mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $SINGULARITY_IMAGEDIR

singularity exec --bind $(pwd):$(pwd) $SINGULARITY_IMAGEDIR/ddh.sif Rscript "$@"
