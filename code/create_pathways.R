library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(vroom)

#rm(list = ls())

go_bp_url <- "https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018"
#consider updating from here: http://amigo.geneontology.org/amigo/software_list

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
  filter(count < 50) %>%  #removes 839 (what is max number of genes you want to query?)
  separate(col = "pathway", into = c("pathway", "go"), sep = "\\(GO\\:") %>% 
  separate(col = "go", into = "go", sep = "\\)", extra = "drop")

rm(go_bp_count)
saveRDS(pathways, file = here::here("data", "pathways.Rds"))


###get something back out
