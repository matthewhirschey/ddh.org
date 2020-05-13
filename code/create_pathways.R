library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(vroom)

#rm(list = ls())

#read current release information to set url for download
source(here::here("code", "current_release.R"))

#GO pathways
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

#GO definitions
go_def <- read_delim(go_def_url, delim = "\n", col_names = "X1") %>% 
  separate(X1, into = c("X1", "X2"), sep = ":", extra = "merge") %>% 
  filter(X1 == "id" | X1 == "name" | X1 == "def") %>% 
  mutate(id = case_when(
      X1 == "id" ~ str_extract(.[[2]], "[:digit:]+"),
      TRUE ~ NA_character_
    )
  ) %>% 
  fill(id) %>% #populate ids across all observations, so that spread can work
  separate(X2, into = "X2", sep = "\\[") %>%  #drop anything after "[", could save it if you want PMIDs, etc.
  select(X1, X2, id) %>% 
  filter(X1 != "id") %>% 
  pivot_wider(names_from = X1, values_from = X2) %>% 
  mutate(def = str_remove_all(def, '\\"'))

#join
pathways <- pathways %>% 
  left_join(go_def, by = c("go" = "id"))
#fix empties (only 28!)
pathways <- pathways %>% 
  mutate(def = case_when(
    is.na(def) ~ "No pathway definition", 
    TRUE ~ def), 
  def = str_trim(def, side = "both")) %>% 
  select(-name)

saveRDS(pathways, file = here::here("data", "pathways.Rds"))
