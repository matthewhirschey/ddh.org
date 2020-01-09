#!/usr/bin/env bash
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
  $RSCRIPT_CMD code/generate_depmap_pathways.R --type $PATHWAY_TYPE --num-subset-files $NUM_SUBSET_FILES --idx $IDX &
done

echo "Waiting for" $PATHWAY_TYPE "pathway processes to finish."
for BGPID in `jobs -p`
do
    wait $BGPID
done

echo "Merging" $PATHWAY_TYPE "pathway subsets."
$RSCRIPT_CMD code/merge_depmap_pathways.R --type $PATHWAY_TYPE --num-subset-files $NUM_SUBSET_FILES
