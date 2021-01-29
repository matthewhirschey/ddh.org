#!/usr/bin/env bash
#SBATCH --job-name=ddh-upload
#SBATCH --mem=16G

if [[ -z "$DDS_PROJECT" ]]
then
    echo "Error: DDS_PROJECT environment variable empty" 1>&2
    exit 1
fi

module load ddsclient
ddsclient upload -p $DDS_PROJECT data
