#!/usr/bin/env bash
#SBATCH --mem=32G
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/ddh-%j.out

source config.sh
export SINGULARITY_TMPDIR SINGULARITY_CACHEDIR ENTREZ_KEY DOCKER_IMG
mkdir -p ${SINGULARITY_TMPDIR} ${SINGULARITY_CACHEDIR} ${SINGULARITY_IMAGEDIR}

make container_image
export RSCRIPT_CMD="singularity exec singularity/images/depmap.sif Rscript"
make
