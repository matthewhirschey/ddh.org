#!/usr/bin/env bash
# Creates a negative or positive pathway file.
# This is done by creating subset files in parallel and then merging them together.
# Environment variables:
# PATHWAY_TYPE determines the pathway file to be created:
#  "positive" ->  master_positive.RData
#  "negative" ->  master_negative.RData
# RSCRIPT_CMD specifies RScript to run
# NUM_SUBSET_FILES specifies number of subset files to create
set -e

if [ -z "$RSCRIPT_CMD" ]; then
  echo "Need to set RSCRIPT_CMD"
  exit 1
fi

if [ -z "$PATHWAY_TYPE" ]; then
  echo "Need to set PATHWAY_TYPE"
  exit 1
fi

if [ -z "$NUM_SUBSET_FILES" ]; then
  echo "Need to set NUM_SUBSET_FILES"
  exit 1
fi

echo "Starting" $NUM_SUBSET_FILES "processes for" $PATHWAY_TYPE "pathway subsets."
for IDX in $(seq $NUM_SUBSET_FILES)
do
  # runs Rscript in the background
  $RSCRIPT_CMD code/generate_depmap_pathways.R --type $PATHWAY_TYPE --num-subset-files $NUM_SUBSET_FILES --idx $IDX &
done

echo "Waiting for" $PATHWAY_TYPE "pathway processes to finish."
for BGPID in `jobs -p`
do
    wait $BGPID
done

echo "Merging" $PATHWAY_TYPE "pathway subsets."
$RSCRIPT_CMD code/merge_depmap_pathways.R --type $PATHWAY_TYPE --num-subset-files $NUM_SUBSET_FILES
