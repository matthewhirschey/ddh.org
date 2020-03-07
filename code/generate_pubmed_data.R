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

#load data
gene_summary <- readRDS(here::here("data", "gene_summary.Rds"))

#download data
pubtatorurl <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz"
gene2pubtator <- read_tsv(pubtatorurl, col_names = c("pmid", "type", "concept_id", "mentions", "resource")) %>% 
  select(pmid, concept_id, mentions) %>% 
  filter(concept_id %in% gene_summary$ncbi_gene_id == TRUE)

gene2pubtator$concept_id <- as.numeric(gene2pubtator$concept_id)

#co-occurance of concepts: count which two genes co-occur in a paper, across all papers
pubmed_concept_pairs <- gene2pubtator %>%
  pairwise_count(concept_id, pmid, sort = TRUE) %>% 
  left_join(gene_summary, by = c("item1" = "ncbi_gene_id")) %>% 
  left_join(gene_summary, by = c("item2" = "ncbi_gene_id")) %>% 
  transmute(target_gene = approved_symbol.x, target_gene_pair = approved_symbol.y, n) %>% 
  pivot_wider(names_from = target_gene_pair, values_from = n, values_fill = list(n = 0)) %>% #this step fills zeros for missing gene-gene pairs
  pivot_longer(names_to = "target_gene_pair", values_to = "n")

#save files
saveRDS(pubmed_concept_pairs, file = here::here("data", paste0(release, "_pubmed_concept_pairs.Rds")))

#how long
time_end_pubmed <- Sys.time()