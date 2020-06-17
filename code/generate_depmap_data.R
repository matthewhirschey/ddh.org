#load libraries
library(tidyverse)
library(here)
library(janitor)
library(corrr)
library(purrr)

#rm(list=ls()) 

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

time_begin_data <- Sys.time()

##BROAD
achilles_raw <- read_csv(achilles_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)")) %>% 
  rename(X1 = 1)

#add name cleaning step
gene_summary <- readRDS(file = here::here("data", "gene_summary.Rds"))
source(here::here("code", "fix_names.R"))
achilles <- clean_colnames(achilles_raw)

achilles_long <- achilles %>% 
  pivot_longer(-X1, names_to = "gene", values_to = "dep_score")

#EXPRESSION(BROAD)
expression <- read_tsv(ccle_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))

#repeat cleaning step for expression
expression <- clean_colnames(expression)

expression_join <- read_csv(cclemeta_url, col_names = TRUE) %>% 
  clean_names() %>% 
  rename(X1 = dep_map_id, cell_line = stripped_cell_line_name) %>% 
  select(X1, cell_line, lineage)

#expression for lineage plots
expression_join_lin <- read_csv(cclemeta_url, col_names = TRUE) %>% 
  clean_names() %>% 
  rename(X1 = dep_map_id, cell_line = stripped_cell_line_name) %>% 
  select(X1, cell_line, lineage, lineage_subtype) %>%
  dplyr::mutate_at("lineage", function(str) {
    str <- str_replace_all(str, "\\_", " ")
    str <- str_to_title(str)
    return(str)
  }) %>%
  dplyr::mutate_at("lineage_subtype", function(str) {
    str <- str_replace_all(str, "\\_", " ")
    str <- str_to_title(str)
    return(str)
  }) 


#filter achilles to remove no expression dep scores(special sauce)
expression_long <- expression %>% 
  filter(expression$X1 %in% achilles$X1 == TRUE) %>% #matches cells
  gather("gene", "gene_expression", -X1) %>% 
  arrange(desc(gene_expression))

no_expression <- expression_long %>% 
  filter(gene_expression == 0) %>% 
  unite(X1, gene, col = "match", sep = "-", remove = TRUE) %>% 
  pull(match)

achilles_no0 <- achilles_long %>% 
  unite(X1, gene, col = "match", sep = "-", remove = FALSE) %>% 
  filter(match %in% no_expression == FALSE) %>% 
  select(-match) %>%
  spread(gene, dep_score)

toomanyNAs <- achilles_no0 %>% 
  summarise_all(list(~sum(is.na(.)))) %>% 
  gather(gene, NAs) %>% 
  arrange(desc(NAs)) %>% 
  filter(NAs > na_cutoff) %>% #set in current_release.R
  pull(gene)

achilles <- achilles_no0 %>% 
  select(-one_of(toomanyNAs)) #check to see if achilles has fewer variables than achilles_raw

#clean Achilles correlation matrix
achilles_cor <- achilles %>%
  select(-X1) %>% 
  correlate() #(diagonal = 0) set to 0 so easy to summarize, but should be NA; so added na.rm = TRUE to fun() in EDA

#save files
saveRDS(achilles, file = here::here("data", paste0(release, "_achilles.Rds")))
saveRDS(expression, file = here::here("data", paste0(release, "_expression.Rds")))
saveRDS(expression_join, file = here::here("data", paste0(release, "_expression_join.Rds")))
saveRDS(expression_join_lin, file = here::here("data", paste0(release, "_expression_join_lin.Rds")))
saveRDS(achilles_cor, file = here::here("data", paste0(release, "_achilles_cor.Rds")))

#how long
time_end_data <- Sys.time()
