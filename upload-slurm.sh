#!/usr/bin/env bash
#SBATCH --mem=16G

module load ddsclient/2.4.0-gcb01
ddsclient upload -p ddh-data data
