#!/usr/bin/env bash
#SBATCH --job-name=ddh-zip
#SBATCH --output=logs/ddh-zip-%j.out
# Zips images and deletes their containing folders.

for IMAGEDIR in "data/images"
do
    BASENAME=$(basename $IMAGEDIR)
    ARCHIVENAME=data/${BASENAME}.tar.gz
    echo "Creating archive $ARCHIVENAME"
    tar -czf $ARCHIVENAME $IMAGEDIR

    echo "Deleting directory $IMAGEDIR"
    rm -rf $IMAGEDIR
done
