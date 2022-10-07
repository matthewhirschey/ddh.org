#!/usr/bin/env bash
#SBATCH --job-name=ddh-pdb
#SBATCH --mem=20G
#SBATCH --output=logs/ddh-pdb-%j.out
set -e

source config.sh
DATADIR=cache/protein_pdbs
PYMOL_CONTAINER=docker://pegi3s/pymol:2.3.0
IMAGEMAGICK_CONTAINER="docker://dpokidov/imagemagick:7.1.0-25"
PROTEIN_PDB_BASE_URL=https://ftp.ebi.ac.uk/pub/databases/alphafold/latest

# Remove github container registry environment variable credentials
unset SINGULARITY_DOCKER_USERNAME 
unset SINGULARITY_DOCKER_PASSWORD

# Ensure $DATDIR exists
mkdir -p $DATADIR

cd $DATADIR 

if [ -e "$PROTEIN_PDB_TARFILE" ]
then
    echo "Using cached $PROTEIN_PDB_BASE_URL/$PROTEIN_PDB_TARFILE file."
else
    echo "Downloading pdb archive $PROTEIN_PDB_BASE_URL/$PROTEIN_PDB_TARFILE"
    wget --quiet $PROTEIN_PDB_BASE_URL/$PROTEIN_PDB_TARFILE
fi

echo "Extracting pdb archive"
tar --overwrite -xf $PROTEIN_PDB_TARFILE

# switch back to ddh base directory
cd ../..

#generate gene list with uniprot IDs
srun slurm/rscript.sh code/generate_gene_csv.R

echo "Convert pdb to png with pymol"
srun singularity exec --pwd /ddh --bind $(pwd):/ddh $PYMOL_CONTAINER pymol -c /ddh/code/generate_png_from_pdb.py $DATADIR/symbol_and_uniprot.csv

echo "Convert png to jpeg with imagemagick mogrify"
srun singularity exec --pwd /ddh --bind $(pwd):/ddh $IMAGEMAGICK_CONTAINER mogrify -format jpg $DATADIR/*.png

