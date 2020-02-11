## THIS CODE GENERATES master_posistive or master_negative TABLES
## Based upon positive_subset_*_of_*.Rds or negative_subset_*_of_*.Rds files created by generate_depmap_pathways.R
## You can run it like so (assuming you ran generate_depmap_pathways.R with --num-subset-files 10)
##  Rscript code/merge_depmap_pathways.R --type positive --num-subset-files 10
##  Rscript code/merge_depmap_pathways.R --type negative --num-subset-files 10

library(tidyverse)
library(here)
library(optparse)

#read dpu_* methods for getting subset filenames and parsing command line parameters
source(here::here("code", "depmap_pathways_util.R"))

master_positive_filename <- "master_positive.Rds"
master_negative_filename <- "master_negative.Rds"

merge_pathways_rds_files <- function(subset_filepaths) {
  merged_data <- tibble(
    fav_gene = character(),
    data = list()
  )
  for (filepath in subset_filepaths) {
    partial_data <- readRDS(filepath)
    merged_data <- merged_data %>%
      bind_rows(partial_data)
  }
  merged_data
}

save_master_positive <- function(subset_filepaths) {
  master_positive <- merge_pathways_rds_files(subset_filepaths)
  message("Saving master_positive to ", master_positive_filename)
  saveRDS(master_positive, file=here::here("data", master_positive_filename))
}

save_master_negative <- function(subset_filepaths) {
  master_negative <- merge_pathways_rds_files(subset_filepaths)
  message("Saving master_negative to ", master_negative_filename)
  saveRDS(master_negative, file=here::here("data", master_negative_filename))
}

main <- function() {
  opt <- dpu_parse_command_line(include_idx=FALSE)
  subset_filepaths <- dpu_get_all_pathways_subset_filepaths(opt$pathways_type, opt$num_subset_files)
  if (opt$pathways_type == dpu_pathways_positive_type) {
    save_master_positive(subset_filepaths)
  } else {
    save_master_negative(subset_filepaths)  
  }
}

if(!interactive()) {
  main()
}
