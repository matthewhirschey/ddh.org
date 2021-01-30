#load libraries
library(tidyverse)
library(here)
library(janitor)
#library(corrr)
#library(moderndive)
#library(purrr)
library(XML)

#rm(list=ls()) 
time_begin_data <- Sys.time()

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#LOAD data 
temp <- tempfile()
download.file(metabolites_url, temp)
hmdb <- xmlToDataFrame(unzip(temp)) %>% #extracting XML took a long time for a 4GB XML file!
  clean_names()
unlink(temp)

hmdb_names <- 
  hmdb %>% 
  select(accession, name, synonyms, pubchem_compound_id)

#save files
saveRDS(hmdb, file = here::here("data", paste0(release, "_hmdb.Rds")))
saveRDS(hmdb_names, file = here::here("data", paste0(release, "_hmdb_names.Rds")))

#how long
time_end_data <- Sys.time()