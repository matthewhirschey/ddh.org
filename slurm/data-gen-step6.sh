#!/usr/bin/env bash
#SBATCH --job-name=ddh-6
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --output=logs/ddh-6-%j.out

set -e

source config.sh


if [ "$GEN_STRUCTURES" == "Y" ]
then
    # Generate images based on pdb files
    slurm/data-gen-pdb.sh

    # Create cache/protein_pdbs.tar.gz for jpg images
    cd cache/protein_pdbs/ \
       && tar -cvf ../protein_pdbs.tar.gz *.jpg \
       && cd ../..

    # Upload cache/protein_pdbs.tar.gz to DDS
    module load ddsclient
    python slurm/dds-upload-file.py $DDS_PROJECT cache/protein_pdbs.tar.gz
else
    # Download cache/protein_pdbs.tar.gz from DDS (if necessary)
    mkdir -p cache
    module load ddsclient
    python slurm/dds-download-file.py $DDS_PROJECT protein_pdbs.tar.gz cache/protein_pdbs.tar.gz
fi

# Download barcodes
python slurm/dds-download-file.py $DDS_PROJECT barcodes.zip cache/barcodes.zip

# Download cell images
python slurm/dds-download-file.py $DDS_PROJECT cell_images.zip cache/cell_images.zip
