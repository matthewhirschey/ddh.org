#~50' to run
#load libraries
library(tidyverse)
library(here)
library(janitor)
library(tidytext)
library(widyr)

#rm(list=ls()) 

time_begin_pubmed <- Sys.time()

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#Load data 
gene_summary <- readRDS(file = here::here("data", "gene_summary.Rds"))
achilles_cor <- readRDS(file = here::here("data", paste0(release, "_achilles_cor.Rds")))
r <- "rowname" #need to drop "rowname"
full <- (names(achilles_cor))[!(names(achilles_cor)) %in% r] #f[!f %in% r]

#download data
#pubtator_url defined in current_release.R
gene2pubtator <- read_tsv(pubtator_url, col_names = c("pmid", "type", "concept_id", "mentions", "resource")) %>% 
  select(pmid, concept_id, mentions) %>% 
  filter(concept_id %in% gene_summary$ncbi_gene_id == TRUE)

gene2pubtator$concept_id <- as.numeric(gene2pubtator$concept_id)

#co-occurrence of concepts: count which two genes co-occur in a paper, across all papers
pubmed_concept_pairs <- gene2pubtator %>%
  pairwise_count(concept_id, pmid, sort = TRUE) %>% 
  left_join(gene_summary, by = c("item1" = "ncbi_gene_id")) %>% 
  left_join(gene_summary, by = c("item2" = "ncbi_gene_id")) %>% 
  transmute(target_gene = approved_symbol.x, target_gene_pair = approved_symbol.y, n)

#nest for faster searching
pubmed_concept_pairs <- pubmed_concept_pairs %>% 
  nest(nested = c(target_gene_pair, n))

#prevent join errors in generate_depmap_tables.R by adding missing genes from names(achilles_cor) to this dataset
missing <- full[which(!full %in% pubmed_concept_pairs$target_gene)]

build_missing_df <- function(missing_gene_vec) {
  missing_dataframe <- tibble(target_gene = character(), 
                              nested = list())
  for (i in missing_gene_vec) {
    mt_df <- tibble(target_gene = as.character(i), 
                    target_gene_pair = as.character(NA), 
                    n = as.numeric(NA))
    
    mt_nest <- mt_df %>% 
      nest(nested = c(target_gene_pair, n))
    
    missing_dataframe <- missing_dataframe %>% 
      bind_rows(mt_nest)
  }
  return(missing_dataframe)
}
missing_df <- build_missing_df(missing)

pubmed_concept_pairs <- pubmed_concept_pairs %>% 
  bind_rows(missing_df)

#save files
saveRDS(pubmed_concept_pairs, file = here::here("data", paste0(release, "_pubmed_concept_pairs.Rds")))

#how long
time_end_pubmed <- Sys.time()
