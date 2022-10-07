#load libraries
library(tidyverse)
library(here)
library(rvest)

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

expression_meta <- readRDS(file=here::here("data", paste0(release, "_expression_meta.Rds")))

#make key
cellosaurus_key <- 
  read_tsv(cellosaurus_url, col_names = FALSE, skip = 34, n_max = 17) %>% 
  tidyr::separate(X1, c("key", "value"), sep = 3)  %>% 
  tidyr::separate(col = value, into = c("value"), sep = "Once|Optional", extra = "drop") %>% 
  dplyr::mutate(key = stringr::str_trim(key, side = "both"), 
                value = stringr::str_trim(value, side = "both"))

#make working df
cellosaurus_raw <- 
  read_tsv(cellosaurus_url, col_names = FALSE, skip = 54) #skip lines which contain data file dictionary

cellosaurus <- 
  cellosaurus_raw %>% 
  tidyr::separate(X1, c("key", "value"), sep = 3) %>% 
  dplyr::mutate(key = stringr::str_trim(key, side = "both"), 
                value = stringr::str_trim(value, side = "both"), 
                id = if_else(grepl("ID", key), value, NA_character_) #must call NA_char so that fill fxn works
  ) %>% 
  tidyr::fill(id) %>% #need fill fxn to populate ids across all observations, so that spread can work
  dplyr::group_by(id, key) %>% 
  dplyr::mutate(values = paste0(value, collapse = "; ")) %>%
  dplyr::select(id, key, values) %>% 
  dplyr::distinct(values, .keep_all = TRUE) %>% 
  tidyr::pivot_wider(names_from = key, values_from = values) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-id)

#add ATCC IDs
cellosaurus <-
  cellosaurus %>%
  dplyr::mutate(ATCC = str_extract(DR, "ATCC; [:graph:]*;")) %>%
  tidyr::separate(ATCC, into = c("tmp", "ATCC"), sep = ";") %>%
  dplyr::mutate(ATCC = stringr::str_trim(ATCC, side = "both")) %>%
  dplyr::select(-tmp)

cellosaurus <-
  cellosaurus %>% 
  dplyr::filter(AC %in% expression_meta$rrid) %>%  #filter 126K cell lines down to only those in CCLE
  dplyr::left_join(expression_meta, by = c("AC" = "rrid")) %>% 
  dplyr::select(name = cell_line, -X1, everything())

saveRDS(cellosaurus_key, file = here::here("data", paste0(release, "_cellosaurus_key.Rds")))
saveRDS(cellosaurus, file = here::here("data", paste0(release, "_cellosaurus.Rds")))
