#LOAD Libs-----
library(tidyverse)
library(furrr)

#LOAD DATA-----
#read current release information
source(here::here("code", "current_release.R"))

#read data from generate_depmap_data.R
achilles <- readRDS(file=here::here(app_data_dir, paste0(release, "_achilles.Rds")))
expression_meta <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression_meta.Rds")))

lineage_counter <- function(expression_data = expression_meta, achilles_data = achilles, gene_symbol) {
  achilles_ids <- 
    achilles_data %>% 
    select(X1, any_of(gene_symbol)) %>% 
    drop_na(gene_symbol) %>% 
    pull(X1)
  
  num <-
    expression_data %>% 
    filter(X1 %in% achilles_ids) %>% 
    distinct(lineage) %>% #lineage_subtype
    nrow(.)
  #nrow(expression_data %>% filter(X1 %in% achilles_ids) %>% count(lineage, sort = TRUE))
  return(num)
}

#test
# lineage_counter(gene_symbol = "RDX")

#map
lineage_summary <- 
  future_map_dfr(colnames(achilles), function(.x) {
  return(data.frame(gene_name = .x, 
                    lineages = lineage_counter(gene_symbol = .x)))
})

lineage_summary %>% 
  count(lineages, sort = TRUE)


