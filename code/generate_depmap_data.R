#load libraries
library(tidyverse)
library(here)
library(janitor)
library(corrr)
library(purrr)

#rm(list=ls()) 

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

start_time <- Sys.time()

##BROAD
achilles <- read_csv(achilles_url, col_names = TRUE) %>%
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))
save(achilles, file = here::here("data", paste0(release, "_achilles.RData")))

achilles_long <- achilles %>% 
  pivot_longer(-X1, names_to = "gene", values_to = "dep_score")

#EXPRESSION(BROAD)
expression <- read_csv(ccle_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))
save(expression, file = here::here("data", paste0(release, "_expression.RData")))

expression_id <- read_csv(cclemeta_url, col_names = TRUE) %>% 
  clean_names
save(expression_id, file = here::here("data", paste0(release, "_expression_id.RData")))

expression_join <- expression_id %>% 
  rename(X1 = dep_map_id) %>% 
  select(X1, stripped_cell_line_name, lineage)
save(expression_join, file = here::here("data", paste0(release, "_expression_join.RData")))

expression_long <- expression %>% 
  filter(expression$X1 %in% achilles$X1 == TRUE) %>% #matches cells
  gather("gene", "gene_expression", -X1) %>% 
  arrange(desc(gene_expression))

no_expression <- expression_long %>% 
  filter(gene_expression == 0) %>% 
  unite(X1, gene, col = "match", sep = "-", remove = TRUE) %>% 
  pull(match) %>% 
  as_tibble()

achilles_no0 <- achilles_long %>% 
  unite(X1, gene, col = "match", sep = "-", remove = FALSE) %>% 
  filter(match %in% no_expression == FALSE) %>% 
  select(-match) %>%
  spread(gene, dep_score) %>% 
  as_tibble()

na_cutoff <- na_cutoff

toomanyNAs <- achilles_no0 %>% 
  summarise_all(list(~sum(is.na(.)))) %>% 
  gather(gene, NAs) %>% 
  arrange(desc(NAs)) %>% 
  filter(NAs > na_cutoff) %>% 
  pull(gene)

achilles_clean <- achilles_no0 %>% 
  select(-one_of(toomanyNAs))

#ACHILLES CORRELATION MATRIX
achilles_cor <- achilles_clean %>% #originally 'achilles'
  select(-X1) %>% 
  correlate() #(diagonal = 0) set to 0 so easy to summarize, but should be NA; so added na.rm = TRUE to fun() in EDA
save(achilles_cor, file = here::here("data", paste0(release, "_achilles_cor.RData")))

#how long
end_time <- Sys.time()
