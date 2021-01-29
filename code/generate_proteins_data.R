#load libraries
library(tidyverse)
library(here)

#rm(list=ls()) 

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#get proteins.Rds from zenodo
tmp_proteins <- tempfile()
download.file(proteins_url, tmp_proteins)
proteins <- readRDS(tmp_proteins)

#save files
saveRDS(proteins, file = here::here("data", paste0(release, "_proteins.Rds")))

