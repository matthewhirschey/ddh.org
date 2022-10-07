#!/usr/bin/env bash
#SBATCH --job-name=ddh-7
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=96G
#SBATCH --output=logs/ddh-7-%j.out

source config.sh
ulimit -s 16384

# Download image.tar.gz from DDS (if necessary)
module load ddsclient
python slurm/dds-download-file.py $DDS_PROJECT images.tar.gz data/images.tar.gz

# Run image generation
DDH_STEP=7 slurm/rscript.sh code/generate_images.R
