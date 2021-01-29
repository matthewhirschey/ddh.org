## THIS CODE GENERATES master_top AND master_bottom TABLES
#~60'

#load libraries
library(tidyverse)
library(here)
library(corrr)

#methods
methods = FALSE

#rm(list=ls()) 
time_begin_tables <- Sys.time()

#read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

#LOAD data 
gene_summary <- readRDS(file = here::here("data", paste0(release, "_gene_summary.Rds")))
achilles_cor <- readRDS(file = here::here("data", paste0(release, "_achilles_cor.Rds")))
achilles_lower <- readRDS(file = here::here("data", paste0(release, "_achilles_lower.Rds")))
achilles_upper <- readRDS(file = here::here("data", paste0(release, "_achilles_upper.Rds")))
mean_virtual_achilles <- readRDS(file = here::here("data", paste0(release, "_mean_virtual_achilles.Rds")))
sd_virtual_achilles <- readRDS(file = here::here("data", paste0(release, "_sd_virtual_achilles.Rds")))
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

#table funs
get_concept_table <- function(concept_table = pubmed_concept_pairs, query_gene) {
  concept_tmp <- concept_table %>% 
    filter(target_gene %in% query_gene) %>% 
    unnest(nested)
  return(concept_tmp)
}

make_dep_table <- function(dep_table = achilles_cor, 
                              summary_table = gene_summary, 
                              query_gene, 
                              upper = achilles_upper, 
                              lower = achilles_lower, 
                              top = TRUE) {
  dep <- 
    dep_table %>% 
    dplyr::select(1, query_gene) %>% 
    {if (top == TRUE) filter(., .[[2]] > upper) else filter(., .[[2]] < lower)} %>% #mean +/- 3sd
    arrange(desc(.[[2]])) %>% #use column index
    left_join(gene_summary, by = c("rowname" = "approved_symbol")) %>% 
    rename(rowname = 1) %>% 
    select(1:3) 
  return(dep)
}

#define list
#sample <- sample(names(achilles_cor), size = 100) #comment this out
r <- "rowname" #need to drop "rowname"
full <- (names(achilles_cor))[!(names(achilles_cor)) %in% r] #f[!f %in% r]
gene_group <- full #(~60' on a laptop); change to sample for testing, or methods

#methods
make_graph_group <- function(graph_gene) {
  top10 <- make_dep_table(query_gene = graph_gene, top = TRUE) %>% 
    top_n(10, wt = .[[2]]) %>% 
    pull(var = 1)
  bottom10 <- make_dep_table(query_gene = graph_gene, top = FALSE) %>% 
    top_n(-10, wt = .[[2]]) %>% 
    pull(var = 1)
  graph_list <- c(graph_gene, top10, bottom10)
  return(graph_list)
}

if(methods == TRUE) {
  methods_gene_query <- make_graph_group(graph_gene = "TP53")
  methods_gene_group <- methods_gene_query
  for (i in methods_gene_query) {
    methods_list <- make_graph_group(graph_gene = i)
    methods_gene_group <- c(methods_gene_group, methods_list)
  }
  gene_group <- unique(methods_gene_group) #this overwrites gene_group with the methods subset
}


#make master_tables
for (fav_gene in gene_group) {
  concept_tmp <- get_concept_table(query_gene = fav_gene)
  
  message(" Dep tables for ", fav_gene)
  dep_top <- make_dep_table(query_gene = fav_gene, top = TRUE)
  dep_top <- 
    dep_top %>% 
    left_join(concept_tmp, by = c("rowname" = "target_gene_pair")) %>% 
    rename(gene = rowname, 
           name = approved_name, 
           r2 = fav_gene,
           concept_count = n) %>% 
    mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles)/sd_virtual_achilles, 1), 
           concept_count = replace_na(concept_count, 0), 
           concept_index = round((concept_count/max(concept_count))*100), 0) %>% 
    select(gene, name, z_score, r2, concept_count, concept_index)
              
    top_table <- dep_top %>% 
      mutate(fav_gene = fav_gene) %>% 
      group_by(fav_gene) %>% 
      nest()
    
    master_top_table <- master_top_table %>% 
      bind_rows(top_table)

  dep_bottom <- make_dep_table(query_gene = fav_gene, top = FALSE)
    
  dep_bottom <-
    dep_bottom %>% 
    left_join(concept_tmp, by = c("rowname" = "target_gene_pair")) %>% 
    rename(gene = rowname, 
           name = approved_name, 
           r2 = fav_gene,
           concept_count = n) %>% 
    mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles)/sd_virtual_achilles, 1), 
           concept_count = replace_na(concept_count, 0), 
           concept_index = round((concept_count/max(concept_count))*100), 0) %>% 
    select(gene, name, z_score, r2, concept_count, concept_index)
  
  bottom_table <- dep_bottom %>% 
    mutate(fav_gene = fav_gene) %>% 
    group_by(fav_gene) %>% 
    nest()
  
  master_bottom_table <- master_bottom_table %>% 
    bind_rows(bottom_table)
}

#save
saveRDS(master_top_table, file=here::here("data", paste0(release, "_master_top_table.Rds")))
saveRDS(master_bottom_table, file=here::here("data", paste0(release, "_master_bottom_table.Rds")))

#make surprise gene list
find_good_candidate <- function(gene_symbol) {
  #this gets the top 10 correlation values
  top_10 <- master_top_table %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(desc(r2)) %>% 
    dplyr::slice(1:10)
  #this looks for 'positive controls' in the top 10...smoking guns
  above <- top_10 %>% 
    dplyr::filter(concept_index > 90) %>% 
    pull(fav_gene)
  #this looks for genes within the top 10 that are under-studied
  below <- top_10 %>% 
    dplyr::filter(concept_index < 10) %>% 
    pull(fav_gene)
  #if both are true, then it's a good candidate for further study
  if(length(above) > 0 && length(below) > 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#TESTING
#find_good_candidate("SSNA1")
#genes <- c("TP53", "TP53BP1")
#map_lgl(genes, ~ find_good_candidate(.))

surprise_genes <- master_top_table %>% 
  mutate(good = purrr::map_lgl(fav_gene, ~ find_good_candidate(.))) %>% 
  dplyr::filter(good == TRUE) %>% 
  pull(fav_gene)

saveRDS(surprise_genes, here::here("data", paste0(release, "_surprise_genes.Rds")))

#Censor
num_genes <- nrow(master_top_table)

genes <- character(num_genes)
num_sim <- numeric(num_genes)

for (i in seq_along(genes)) {
  genes[i] <- master_top_table$fav_gene[i]
  num_sim[i] <- nrow(master_top_table[[2]][[i]])
}

censor_genes <- tibble(genes, num_sim)

#save
saveRDS(censor_genes, here::here("data", paste0(release, "_censor_genes.Rds")))

#how long
time_end_tables <- Sys.time()
