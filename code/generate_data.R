## Generates all data needed for running the shiny app one item at a time.
## To generate data faster use a slurm cluster see the submit-slurm-jobs.sh script
## Uses DDH_STEP environment variable to determine which step to run.
## Uses DDH_IDX environment variable when DDH_STEP=4.pos or DDH_STEP=4.neg

library(here)
source(here::here("code", "current_release.R"))
source(here::here("code", "install_libraries.R"))
# check DDH_PRIVATE to determine if we are generating public vs private data
private = Sys.getenv("DDH_PRIVATE") == "Y"

generate_data_step1 <- function() {
  # Requires 64G of memory
  message("DDH: Running step 1 - Gene summary and other files.")
  #SUMMARY
  source(here::here("code", "generate_gene_summary.R"))
  generate_gene_summary(gene_names_url, entrez_key, here::here("data", paste0(release, "_", gene_summary_output_filename)), gene_symbol = NULL, regenerate = FALSE)
  message("generate gene summary complete")
  #PATHWAYS
  source(here::here("code", "generate_pathways.R"))
  message("generate pathways complete")
  #GENE LOCATION & CHROMOSOME
  source(here::here("code", "generate_gene_location.R"))
  generate_gene_location(regenerate = FALSE)
  message("generate gene location complete")
  #THRESHOLD
  source(here::here("code", "generate_threshold.R")) #run this first to generate the na_cutoff
  generate_threshold()
  message("generate threshold complete")
  #DEPMAP
  source(here::here("code", "generate_depmap_data.R")) #script to generate ddh correlation matrix
  # source(here::here("code", "generate_GLS.R")) #script to generate ddh GLS pvalues
  message("generate depmap data complete")
  #SUBCELL
  source(here::here("code", "generate_subcell_data.R")) #script to generate subcell data
  message("generate subcell data complete")
  #TISSUE
  source(here::here("code", "generate_tissue_data.R")) #script to generate tissue data
  message("generate tissue data complete")
  #PROTEINS
  if (private) {
    source(here::here("code", "generate_proteins_data_private.R")) #script to generate proteins data
  } else {
    source(here::here("code", "generate_proteins_data.R")) #script to generate proteins data
  }
  ##SEQUENCE CLUSTERS
  message("generating protein sequence clusters")
  source(here::here("code", "generate_sequence_clusters.R"))
  message("generating protein clusters enrichment")
  source(here::here("code", "generate_protein_cluster_enrichment.R"))
  message("generate proteins data complete")
  ##PROTEIN DOMAINS
  source(here::here("code", "generate_protein_domains.R"))
  message("generate protein domain data complete")
  ##PDB IDs
  source(here::here("code", "generate_pdb_ids.R"))
  message("generate PDB IDs complete")
  #STATS
  source(here::here("code", "generate_depmap_stats.R")) #script to generate ddh stats
  message("generate stats complete")
  
  message("DDH: Finished step 1.")
}

generate_data_step2 <- function() {
  # Requires 212G of memory
  message("DDH: Running step 2 - drug & pubmed data.")
  #script to generate drug data
  if (private) {
    source(here::here("code", "generate_drug_data_private.R"))
  } else {
    source(here::here("code", "generate_drug_data.R"))
  }
  source(here::here("code", "generate_pubmed_data.R")) #script to generate pubtator relationships
  source(here::here("code", "generate_metabolites_data.R")) #generates metabolite data and collapses
  #CELL LINE (needs drug data)
  if (private) {
    source(here::here("code", "generate_cell_line_data_private.R")) #script to generate cell line data
  } else {
    source(here::here("code", "generate_cell_line_data.R")) 
  }
  message("generate cell line complete")
  #CELL LINE SIMILARITY
  source(here::here("code", "generate_cell_line_similarity.R"))
  message("generate cell line similarity complete")
  
  message("DDH: Finished step 2.")
}

generate_data_step3 <- function() {
  # Requires 32G of memory
  message("DDH: Running step 3 - ddh tables.")
  source(here::here("code", "generate_depmap_tables.R")) #third script to generate ddh tables
  source(here::here("code", "generate_surprise.R"))
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

generate_data_step5 <- function() {
  # Requires 112G of memory
  message("DDH: Running step 5 - starting DDH methods")
    source(here::here("code", "generate_methods.R")) #data should be there, so just need to render
  message("DDH: Finished step 5.")
}

generate_data_step6 <- function() {
  # Requires 96G of memory <- NOT SURE ABOUT THIS;
  message("DDH: Running step 6 - starting structure image generation")
  source(here::here("code", "generate_structures.R"))
  message("DDH: Finished step 6.")
}

generate_data_step7 <- function() {
  # Requires 96G of memory <- NOT SURE ABOUT THIS; UPDATE + slurm/data-gen-step6.sh
  message("DDH: Running step 7 - starting image generation")
  source(here::here("code", "generate_images.R"))
  message("DDH: Finished step 7.")
}

generate_data <- function() {
  message("Generating data for ", release, ".")
  generate_data_step1()
  generate_data_step2()
  generate_data_step3()
  generate_data_step4()
  generate_data_step5()
  generate_data_step6()
  generate_data_step7()
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
  "4.merge"=generate_data_step4_merge_subsets, 
  "5"=generate_data_step5,
  "6"=generate_data_step6, 
  "7"=generate_data_step7
)
ddh_step_to_func[[ddh_step]]()
