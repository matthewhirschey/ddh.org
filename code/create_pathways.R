library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(vroom)

#rm(list = ls())

#read current release information to set url for download
source(here::here("code", "current_release.R"))

go_bp <- vroom(go_bp_url, col_names = FALSE, delim = "\t") %>%  #, col_names = FALSE, .name_repair = "universal") %>% 
  clean_names() %>% 
  select(-x2) %>% 
  unite(col = "gene", x3:x15, sep = "\t") %>% 
  separate_rows(gene, sep = "\t") %>% 
  rename("pathway" = "x1") %>% 
  filter(!is.na(gene), 
         gene != "NA", 
         gene != "") %>% 
  group_by(pathway)

go_bp_count <- go_bp %>%  
  summarize(count = n()) %>% 
  ungroup()

pathways <- go_bp %>% 
  nest() %>% 
  left_join(go_bp_count, by = "pathway") %>% 
  filter(count < 41) %>%  #5103 total; <50 removes 839; <30 removes 1462; <20 removes 2118
  separate(col = "pathway", into = c("pathway", "go"), sep = "\\(GO\\:") %>% 
  separate(col = "go", into = "go", sep = "\\)", extra = "drop") %>% 
  mutate(pathway = str_trim(pathway, side = "right"), 
         pathway = str_to_title(pathway), 
         pathway = str_replace_all(pathway, "Ii", "II"),
         pathway = str_replace_all(pathway, "ii", "II"))

rm(go_bp_count)
saveRDS(pathways, file = here::here("data", "pathways.Rds"))
