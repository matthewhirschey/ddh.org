##This script updates gene_summary to include pubmed data

#load libraries
library(tidyverse)
library(here)
library(janitor)

#rm(list=ls()) 

time_begin_gene_summary<- Sys.time()

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#import gene_summary, and select original columns
gene_summary <- readRDS(here::here("data", "gene_summary_raw.Rds")) %>% 
  select ("approved_symbol", "approved_name", "aka", "ncbi_gene_id", "hgnc_id", "chromosome", "ref_seq_i_ds", "locus_type", "omim_id_supplied_by_omim", "uni_prot_id_supplied_by_uni_prot", "entrez_summary")

#gene2pubmedurl defined in current_release.R
gene2pubmed_raw <- read_tsv(gene2pubmed_url, col_names = TRUE) %>% 
  clean_names()

#add pubmed_count to gene_summary
gene_summary <- gene2pubmed_raw %>% 
  filter(number_tax_id == 9606) %>%  #only the rows corresponding to humans (#tax_id = 9606) 
  group_by(gene_id) %>% 
  count(sort = TRUE) %>% 
  right_join(gene_summary, by = c("gene_id" = "ncbi_gene_id"))  %>% 
  rename("ncbi_gene_id" = "gene_id", 
         "pubmed_count" = "n") %>% 
  mutate(pubmed_count = replace_na(pubmed_count, 0)) %>% 
  ungroup()

#add concept count
pubtatorurl <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz"
gene2pubtator <- read_tsv(pubtatorurl, col_names = c("pmid", "type", "concept_id", "mentions", "resource")) %>% 
  select(pmid, concept_id, mentions) %>% 
  filter(concept_id %in% gene_summary$ncbi_gene_id == TRUE)

gene2pubtator$concept_id <- as.numeric(gene2pubtator$concept_id)

gene_summary <- gene2pubtator %>% 
  group_by(concept_id) %>% 
  count(sort = TRUE) %>% 
  right_join(gene_summary, by = c("concept_id" = "ncbi_gene_id"))  %>% 
  rename("ncbi_gene_id" = "concept_id", 
         "concept_count" = "n") %>% 
  mutate(concept_count = replace_na(concept_count, 0)) %>% 
  ungroup()

gene_summary <- gene_summary %>%
  mutate(pubmed_count_rank = percent_rank(pubmed_count),
         concept_count_rank = percent_rank(concept_count)) %>% 
  select("approved_symbol", "approved_name", "aka", "ncbi_gene_id", "hgnc_id", "chromosome", "ref_seq_i_ds", "locus_type", "omim_id_supplied_by_omim", "uni_prot_id_supplied_by_uni_prot", "entrez_summary", "pubmed_count", "pubmed_count_rank", "concept_count", "concept_count_rank")

#save final file
saveRDS(gene_summary, here::here("data", "gene_summary.Rds"))

#how long
time_end_gene_summary <- Sys.time()
