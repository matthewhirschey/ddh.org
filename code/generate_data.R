## Generates all data needed for running the shiny app one item at a time.
## To generate data faster use a slurm cluster see the submit-slurm-jobs.sh script
## Uses DDH_STEP environment variable to determine which step to run.
## Uses DDH_IDX environment variable when DDH_STEP=4.pos or DDH_STEP=4.neg

library(here)
source(here::here("code", "current_release.R"))

generate_data_step1 <- function() {
  # Requires 64G of memory
  message("DDH: Running step 1 - Gene summary and other files.")
  source(here::here("code", "create_gene_summary.R"))
  source(here::here("code", "create_pathways.R"))
  source(here::here("code", "find_threshold.R")) #run this first to generate the na_cutoff
  source(here::here("code", "generate_depmap_data.R")) #script to generate ddh correlation matrix
  source(here::here("code", "generate_depmap_stats.R")) #script to generate ddh stats
  source(here::here("code", "generate_subcell_data.R")) #script to generate subcell data
  source(here::here("code", "generate_proteins_data.R")) #script to generate proteins data
  message("DDH: Finished step 1.")
}

generate_data_step2 <- function() {
  # Requires 212G of memory
  message("DDH: Running step 2 - pubmed data.")
  source(here::here("code", "generate_drug_data_names.R")) #script to generate proteins data
  source(here::here("code", "generate_pubmed_data.R")) #script to generate pubtator relationships
  message("DDH: Finished step 2.")
}

generate_data_step3 <- function() {
  # Requires 32G of memory
  message("DDH: Running step 3 - ddh tables.")
  source(here::here("code", "generate_depmap_tables.R")) #third script to generate ddh tables
  message("DDH: Finished step 3.")
}

generate_data_step4_positive_subset <- function(idx=strtoi(Sys.getenv("DDH_IDX"))) {
  message("DDH: Running step 4.pos ", idx)
  source(here::here("code", "generate_depmap_pathways.R"))
  if (!file.exists(dpu_get_subset_filepath(dpu_pathways_positive_type, idx, pathways_num_subset_files))) {
    create_pathway_subset_file(dpu_pathways_positive_type, idx, pathways_num_subset_files)
  }
  message("DDH: Finished step 4.pos ", idx)
}

generate_data_step4_negative_subset <- function(idx=strtoi(Sys.getenv("DDH_IDX"))) {
  message("DDH: Running step 4.neg ", idx)
  source(here::here("code", "generate_depmap_pathways.R"))
  if (!file.exists(dpu_get_subset_filepath(dpu_pathways_negative_type, idx, pathways_num_subset_files))) {
    create_pathway_subset_file(dpu_pathways_negative_type, idx, pathways_num_subset_files)
  }
  message("DDH: Finished step 4.neg ", idx)
}

generate_data_step4_merge_subsets <- function() {
  source(here::here("code", "generate_depmap_pathways.R"))
  message("DDH: Running step 4.merge")
  save_master_positive(dpu_get_all_pathways_subset_filepaths(dpu_pathways_positive_type, pathways_num_subset_files))
  save_master_negative(dpu_get_all_pathways_subset_filepaths(dpu_pathways_negative_type, pathways_num_subset_files))
  message("DDH: Deleting pathway subset files.")
  for(i in 1:pathways_num_subset_files) {
    unlink(dpu_get_subset_filepath(dpu_pathways_positive_type, i, pathways_num_subset_files))
    unlink(dpu_get_subset_filepath(dpu_pathways_negative_type, i, pathways_num_subset_files))
  }
  message("DDH: Finished step 4.merge")
}

generate_data_step4 <- function() {
  message("DDH: Running step 4 - pathway data generation.")
  for(idx in 1:pathways_num_subset_files) {
    generate_data_step4_positive_subset(idx)
  }
  for(idx in 1:pathways_num_subset_files) {
    generate_data_step4_negative_subset(idx)
  }
  generate_data_step4_merge_subsets()
  message("DDH: Finished step 4 - pathway data generation.")
}

generate_data <- function() {
  message("Generating data for ", release, ".")
  generate_data_step1()
  generate_data_step2()
  generate_data_step3()
  generate_data_step4()
  message("Done.")
}

# Read DDH_STEP environment variable and run the requested step
# Defaults to running all steps if the environment variable is not set
ddh_step <- Sys.getenv("DDH_STEP", "ALL")
ddh_step_to_func <- list(
  "ALL"=generate_data,
  "1"=generate_data_step1,
  "2"=generate_data_step2,
  "3"=generate_data_step3,
  "4.pos"=generate_data_step4_positive_subset,
  "4.neg"=generate_data_step4_negative_subset,
  "4.merge"=generate_data_step4_merge_subsets
)
ddh_step_to_func[[ddh_step]]()
