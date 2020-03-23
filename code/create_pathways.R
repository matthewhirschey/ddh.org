library(tidyverse)
library(here)
library(readxl)
library(janitor)

go_bp <- read_xls(here::here("data", "GO_Biological_Process_2018.xls"), col_names = FALSE, .name_repair = "universal") %>% 
  clean_names() %>%  
  select(-x2) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = -x1, names_to = "gene", values_to = "gene_name") %>% 
  filter(!is.na(gene_name)) %>% 
  select(-gene) %>% 
  rename("pathway" = "x1") %>% 
  group_by(pathway)

go_bp_count <- go_bp %>%  
  summarize(count = n()) %>% 
  ungroup()

go_bp <- go_bp %>% 
  nest() %>% 
  left_join(go_bp_count, by = "pathway")
  
