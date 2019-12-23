# Utility functions for use with depmap pathway generation
library(here)

dpu_pathways_positive_type <- "positive"
dpu_pathways_negative_type <- "negative"

# Get a path to a subset file within the data directory for either positive or negative pathway data
dpu_get_subset_filepath <- function(pathways_type, file_idx, num_subset_files) {
  # create filename like 'positive_subset_1_of_10' or 'negative_subset_1_of_10'
  subset_filename <- paste0(pathways_type, "_subset", "_", file_idx, "_of_", num_subset_files, ".Rds")
  here::here("data", subset_filename)
}

# Get a path all subset files within the data directory for either positive or negative pathway data
dpu_get_all_pathways_subset_filepaths <- function(pathways_type, num_subset_files) {
  dpu_get_subset_filepath(pathways_type, seq(num_subset_files), num_subset_files)
}

# parse command line for --type, --num-subset-files and optionally --idx
dpu_parse_command_line <- function (include_idx) {
  option_list <- list(
    make_option(c("--type"), type="character", default="positive", dest='pathways_type',
                help="Type of pathways file to create either 'positive' or 'negative'"),
    make_option(c("--num-subset-files"), type="integer", dest='num_subset_files',
                help="Number of subset pathways files")
  )
  if (include_idx) {
    option_list <- c(
      option_list,
      make_option(c("--idx"), type="integer", dest='subset_file_idx',
                  help="Specifies single pathways subset file to create.")
    )
  }
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  if (opt$pathways_type != dpu_pathways_positive_type && opt$pathways_type != dpu_pathways_negative_type) {
    message("The --type argument must be either '", dpu_pathways_positive_type, "' or '", dpu_pathways_positive_type, "'")
    q(status=1)
  }
  if (is.null(opt$num_subset_files)) {
    message("The --num-subset-files flag is required.")
    q(status=1)
  }
  if (include_idx && is.null(opt$subset_file_idx)) {
    message("The --idx is required.")
    q(status=1)
  }
  opt
}