library(here)
source(here::here("code", "current_release.R"))

message("Generating data for ", release, ".")

# read DDH_STEP environment variable to get steps to run or the default which is all steps.
steps_to_run = unlist(strsplit(Sys.getenv("DDH_STEP", "1,2,3,4"), ','))
# read ENTREZ_KEY from environment variable
entrez_key = Sys.getenv("ENTREZ_KEY")

if ("1" %in% steps_to_run) {
  message("DDH: Running step 1.")
  source(here::here("code", "create_gene_summary.R"))
  create_gene_summary(gene_names_url, entrez_key, here::here("data", gene_summary_output_filename))

  source(here::here("code", "find_threshold.R")) #run this first to generate the na_cutoff
  source(here::here("code", "generate_depmap_data.R")) #script to generate ddh correlation matrix
  source(here::here("code", "generate_depmap_stats.R")) #script to generate ddh stats
  message("DDH: Finished step 1.")
}

if ("2" %in% steps_to_run) {
  message("DDH: Running step 2.")
  source(here::here("code", "generate_pubmed_data.R")) #script to generate pubtator relationships
  message("DDH: Finished step 2.")
}

if ("3" %in% steps_to_run) {
  message("DDH: Running step 3.")
  source(here::here("code", "generate_depmap_tables.R")) #third script to generate ddh tables
  message("DDH: Finished step 3.")
}

if ("4" %in% steps_to_run) {
  message("DDH: Running step 4.")
  source(here::here("code", "generate_depmap_pathways.R"))
  num_subset_files <- 10

  library(doParallel)
  cl <- makeCluster(num_subset_files)
  registerDoParallel(cl)
  dpu_pathways_positive_type <- "positive"
  dpu_pathways_negative_type <- "negative"

  foreach(i=1:num_subset_files) %dopar% create_pathway_subset_file(dpu_pathways_positive_type, i, num_subset_files)
  save_master_positive(subset_filepaths)

  foreach(i=1:num_subset_files) %dopar% create_pathway_subset_file(dpu_pathways_negative_type, 2, num_subset_files)
  save_master_negative(subset_filepaths)

  stopCluster(cl)
  message("DDH: Finished step 4.")
}

message("Done.")

