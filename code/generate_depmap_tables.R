## THIS CODE GENERATES master_top AND master_bottom TABLES
#~60'

#load libraries
library(tidyverse)
library(here)
library(corrr)

#rm(list=ls()) 
time_begin_tables <- Sys.time()

#read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

#LOAD data 
gene_summary <- readRDS(file = here::here("data", "gene_summary.Rds"))
achilles_cor <- readRDS(file = here::here("data", paste0(release, "_achilles_cor.Rds")))
achilles_lower <- readRDS(file = here::here("data", "achilles_lower.Rds"))
achilles_upper <- readRDS(file = here::here("data", "achilles_upper.Rds"))
mean_virtual_achilles <- readRDS(file = here::here("data", "mean_virtual_achilles.Rds"))
sd_virtual_achilles <- readRDS(file = here::here("data", "sd_virtual_achilles.Rds"))
pubmed_concept_pairs <- readRDS(file = here::here("data", paste0(release, "_pubmed_concept_pairs.Rds")))

#setup containers
master_top_table <- tibble(
  fav_gene = character(), 
  data = list()
)
master_bottom_table <- tibble(
  fav_gene = character(), 
  data = list()
)

#define list
#sample <- sample(names(achilles_cor), size = 100) #comment this out
r <- "rowname" #need to drop "rowname"
full <- (names(achilles_cor))[!(names(achilles_cor)) %in% r] #f[!f %in% r]

gene_group <- full #(~60' on a laptop); change to sample for testing

#master_tables
for (fav_gene in gene_group) {
  concept_tmp <- pubmed_concept_pairs %>% 
    filter(target_gene == fav_gene) %>% 
    unnest(nested)
  
  message(" Dep tables for ", fav_gene)
  dep_top <- achilles_cor %>% 
    focus(fav_gene) %>% 
    arrange(desc(.[[2]])) %>% #use column index
    filter(.[[2]] > achilles_upper) %>% #mean +/- 3sd
    left_join(gene_summary, by = c("rowname" = "approved_symbol")) %>% 
    select(rowname, approved_name, fav_gene) %>% 
    left_join(concept_tmp, by = c("rowname" = "target_gene_pair")) %>% 
    rename(gene = rowname, 
           name = approved_name, 
           r2 = fav_gene,
           concept_count = n) %>% 
    mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles)/sd_virtual_achilles, 1), 
           concept_count = replace_na(concept_count, 0)) %>% 
    select(gene, name, z_score, r2, concept_count)
              
    top_table <- dep_top %>% 
      mutate(fav_gene = fav_gene) %>% 
      group_by(fav_gene) %>% 
      nest()
    
    master_top_table <- master_top_table %>% 
      bind_rows(top_table)

  dep_bottom <- achilles_cor %>% 
    focus(fav_gene) %>% 
    arrange(desc(.[[2]])) %>% #use column index
    filter(.[[2]] < achilles_lower) %>% #mean +/- 3sd
    left_join(gene_summary, by = c("rowname" = "approved_symbol")) %>% 
    select(rowname, approved_name, fav_gene) %>% 
    left_join(concept_tmp, by = c("rowname" = "target_gene_pair")) %>% 
    rename(gene = rowname, 
           name = approved_name, 
           r2 = fav_gene,
           concept_count = n) %>% 
    mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles)/sd_virtual_achilles, 1), 
           concept_count = replace_na(concept_count, 0)) %>% 
    select(gene, name, z_score, r2, concept_count)
  
  bottom_table <- dep_bottom %>% 
    mutate(fav_gene = fav_gene) %>% 
    group_by(fav_gene) %>% 
    nest()
  
  master_bottom_table <- master_bottom_table %>% 
    bind_rows(bottom_table)
}

#save
saveRDS(master_top_table, file=here::here("data", "master_top_table.Rds"))
saveRDS(master_bottom_table, file=here::here("data", "master_bottom_table.Rds"))

#how long
time_end_tables <- Sys.time()
