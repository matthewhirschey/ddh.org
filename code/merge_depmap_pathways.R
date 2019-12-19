library(tidyverse)
library(here)
library(optparse)

merge_rds_files <- function(input_filename_prefix, group_size) {
  merged_data <- tibble(
    fav_gene = character(),
    data = list()
  )
  for (group_idx in seq(group_size)) {
    output_filename_suffix <- paste0("_", group_idx, "_of_", group_size)
    input_filename <- paste0(input_filename_prefix, output_filename_suffix, ".Rds")
    partial_data <- readRDS(here::here("data", input_filename))
    merged_data <- merged_data %>%
      bind_rows(partial_data)
  }
  merged_data
}

option_list <- list(
  make_option(c("--groups"), type="integer", default=1,
              help="Number of groups to split genes into")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
group_size <- opt$groups

master_positive_filename <- "master_positive.RData"
master_positive <- merge_rds_files("master_positive", group_size)
print(paste0("Saving master_positive to ", master_positive_filename))
saveRDS(master_positive, file=here::here("data", master_positive_filename))

master_negative_filename <- "master_negative.RData"
master_negative <- merge_rds_files("master_negative", group_size)
print(paste0("Saving master_negative to ", master_negative_filename))
saveRDS(master_negative, file=here::here("data", master_negative_filename))
