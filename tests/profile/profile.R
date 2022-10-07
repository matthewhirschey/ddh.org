library(tictoc)
library(pryr)
library(shiny)
source(here::here("code", "install_libraries.R"))

#LOAD DATA-----
#read params
source(here::here("code", "app_params.R"))
source(here::here("code", "generate_colors.R"))

# replace readRDS with a version that tracks how long each file takes to load
sv_readRDS <- readRDS

readRDS <- function(file) {
  tictoc::tic(basename(file))
  res <- sv_readRDS(file)
  tictoc::toc(log = TRUE, quiet = TRUE)
  res
}

message("Loading data files...")

tictoc::tic.clearlog()
source(here::here("code", "app_data.R"))
app_data_log_list <- tictoc::tic.log(format=FALSE)

log_df <- data.frame(
  elapsed=vapply(app_data_log_list, function(x) x$toc - x$tic, 0),
  msg=vapply(app_data_log_list, function(x) x$msg, '')
)
sorted_log_df <- dplyr::arrange(log_df, desc(elapsed))
app_data_top_10 <- head(sorted_log_df, 10)

tictoc::tic.clearlog()

testShinyModule <- function(moduleName, ...) {
  print(moduleName)
  tictoc::tic(moduleName)
  tryCatch({
    testServer(get(moduleName), ...)
  }, error=function(e){
    print("ERROR")
    print(e)
  })
  tictoc::toc(log = TRUE, quiet = TRUE)
  invisible()
}

source(here::here("tests/profile", "profile_args.R"))

print("Profiling shiny_text.R")
source(here::here("tests/profile", "profile_shiny_text.R"))

print("Profiling shiny_plots.R")
source(here::here("tests/profile", "profile_shiny_plots.R"))

print("Profiling shiny_graphs.R")
source(here::here("tests/profile", "profile_shiny_graphs.R"))

print("Profiling shiny_tables.R")
source(here::here("tests/profile", "profile_shiny_tables.R"))

print("Profiling shiny_cards.R")
source(here::here("tests/profile", "profile_shiny_cards.R"))

print("Profiling shiny_structures.R")
source(here::here("tests/profile", "profile_shiny_structures.R"))

print("Module Results:")
writeLines(unlist(tictoc::tic.log()))

message("")
message("Top 50 slowest modules:")
log_list <- tictoc::tic.log(format=FALSE)
log_df <- data.frame(
  elapsed=vapply(log_list, function(x) x$toc - x$tic, 0),
  msg=vapply(log_list, function(x) x$msg, '')
)
sorted_log_df <- dplyr::arrange(log_df, desc(elapsed))
print(head(sorted_log_df, 50))

message("")
message("Top 10 slowest loading files:")
print(app_data_top_10)

message("Memory used:")
print(mem_used())
