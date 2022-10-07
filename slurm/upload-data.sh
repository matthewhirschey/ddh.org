#!/usr/bin/env bash
#SBATCH --job-name=ddh-upload
#SBATCH --mem=16G
#SBATCH --output=logs/ddh-upload-data-%j.out

if [[ -z "$DDS_PROJECT" ]]
then
    echo "Error: DDS_PROJECT environment variable empty" 1>&2
    exit 1
fi
set -e

module load ddsclient
ddsclient upload -p $DDS_PROJECT data

if [ "$DELETE_DATA_AFTER_UPLOAD" = "Y" ]
then
   rm -rf data
fi
