#!/usr/bin/env bash
#SBATCH --mem=32G


# 16G sufficient for depmap data
# asking 32 for stats
export SINGULARITY_TMPDIR=/data/itlab/singularity_tmp
export SINGULARITY_CACHEDIR=/data/itlab/singularity_cache

make all
