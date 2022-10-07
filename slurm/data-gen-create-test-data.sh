#!/usr/bin/env bash
#SBATCH --job-name=ddh-testdata
#SBATCH --output=logs/ddh-testdata-%j.out
#SBATCH --mem=10G
set -e

source config.sh

echo "Deleting old test data directory"
rm tests/data/*.Rds

echo "Running create_test_data.R"
slurm/rscript.sh code/create_test_data.R

echo "Done"
