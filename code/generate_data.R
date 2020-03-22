## Generates all data needed for running the shiny app
## Environment variable DDH_STEP can be used to control which steps are run.
## The step names are "gens", "pubm", "tbls", "path". They need to be run in this order.
## By default all steps are run in the correct order. These steps are broken up based on 
## Their different resource requirments.

library(here)
source(here::here("code", "current_release.R"))

message("Generating data for ", release, ".")

# read DDH_STEP environment variable to get steps to run or the default which is all steps.
steps_to_run = unlist(strsplit(Sys.getenv("DDH_STEP", "gens,pubm,tbls,path"), ','))
# read ENTREZ_KEY from environment variable
entrez_key = Sys.getenv("ENTREZ_KEY")

if ("gens" %in% steps_to_run) {
  # Requires around 64G of memory
  message("DDH: Running step 1 - Gene summary and other files.")
  source(here::here("code", "create_gene_summary.R"))
  source(here::here("code", "find_threshold.R")) #run this first to generate the na_cutoff
  source(here::here("code", "generate_depmap_data.R")) #script to generate ddh correlation matrix
  source(here::here("code", "generate_depmap_stats.R")) #script to generate ddh stats
  message("DDH: Finished step 1.")
}

if ("pubm" %in% steps_to_run) {
  # Requires around 212G of memory
  message("DDH: Running step 2 - pubmed data.")
  source(here::here("code", "generate_pubmed_data.R")) #script to generate pubtator relationships
  message("DDH: Finished step 2.")
}

if ("tbls" %in% steps_to_run) {
  # Requires 32G of memory
  message("DDH: Running step 3 - ddh tables.")
  source(here::here("code", "generate_depmap_tables.R")) #third script to generate ddh tables
  message("DDH: Finished step 3.")
}

if ("path" %in% steps_to_run) {
  # Requires 32G of memory and 10 CPUs
  message("DDH: Running step 4 - pathway data generation.")
  source(here::here("code", "generate_depmap_pathways.R"))
  num_subset_files <- 10

  library(doParallel)
  cl <- makeCluster(num_subset_files)
  registerDoParallel(cl)
  dpu_pathways_positive_type <- "positive"
  dpu_pathways_negative_type <- "negative"

  message("DDH: Creating pathway postive subset files.")
  foreach(i=1:num_subset_files) %dopar% {
      source(here::here("code", "generate_depmap_pathways.R"))
      create_pathway_subset_file(dpu_pathways_positive_type, i, num_subset_files)
  }
  message("DDH: Save master postive file.")
  save_master_positive(dpu_get_all_pathways_subset_filepaths(dpu_pathways_positive_type, num_subset_files))

  message("DDH: Creating pathway negative subset files.")
  foreach(i=1:num_subset_files) %dopar% {
      source(here::here("code", "generate_depmap_pathways.R"))
      create_pathway_subset_file(dpu_pathways_negative_type, i, num_subset_files)
  }
  message("DDH: Save master negative file.")
  save_master_negative(dpu_get_all_pathways_subset_filepaths(dpu_pathways_negative_type, num_subset_files))

  stopCluster(cl)
  message("DDH: Finished step 4.")
}

message("Done.")

