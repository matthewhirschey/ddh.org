#!/usr/bin/env bash
#SBATCH --mem=16G
#SBATCH --output=logs/upload-%j.out

module load ddsclient/2.4.0-gcb01

source config.sh
cd ${DDH_BASE}
ddsclient upload -p ddh-data data
