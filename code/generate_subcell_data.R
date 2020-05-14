library(tidyverse)
library(here)
library(janitor)
library(gganatogram)

source(here::here("code", "current_release.R"))

#import
zip_file <- tempfile(fileext = ".zip")
download.file(subcell_url, zip_file, mode = "wb")
subcell <- read_tsv(zip_file) %>% 
  clean_names()

#cleaning
subcell <- subcell %>% 
  separate_rows(main_location, sep = "\\;")

subcell$organ <- str_to_lower(subcell$main_location)
subcell$organ <- str_replace(subcell$organ, " ", "_")

#join
data(cell_key)

subcell <- subcell %>% 
  full_join(cell_key[['cell']], by = "organ")

#save file
saveRDS(subcell, file = here::here("data", paste0(release, "_subcell.Rds")))
